/*
 *                             The MIT License
 *
 * Copyright (c) 2024 by Albert Jimenez-Blanco
 *
 * This file is part of #################### Theseus Library ####################.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */


#include <string_view>
#include <chrono>
#include <atomic>
#include <queue>
#include <algorithm>
#include <cstdint>
#include "theseus_aligner_impl.h"

// Cumulative timing counters for align_from profiling.
std::atomic<int64_t> g_subgraph_ns{0};
std::atomic<int64_t> g_new_alignment_ns{0};
std::atomic<int64_t> g_dp_loop_ns{0};
std::atomic<int64_t> g_backtrace_ns{0};
std::atomic<int64_t> g_poa_update_ns{0};
std::atomic<int64_t> g_align_calls{0};
std::atomic<int64_t> g_total_score{0};

namespace theseus {

TheseusAlignerImpl::TheseusAlignerImpl(const Penalties &penalties,
                                       const Heuristics &heuristics,
                                       Graph &&graph,
                                       int initial_weight,
                                       NodeId sink_node) :    _penalties(penalties),
                                                              _internal_penalties(penalties),
                                                              _heuristics(heuristics),
                                                              _graph(std::move(graph)),
                                                              _seq("", false) {
    _source_node = 0;
    _sink_node = sink_node;
    // TODO: Gap-linear and dual affine-gap.
    const auto n_scores = std::max({_internal_penalties.gapo() +_internal_penalties.gape(),
                                  _internal_penalties.gapo() +_internal_penalties.gape(),
                                  _internal_penalties.mism()}) + 1;

    // Always MSA for pericles
    _poa_graph = std::make_unique<POAGraph>();
    if (_sink_node == 2) {
        _poa_graph->create_initial_graph(_graph, initial_weight);
    }
    else {
        _poa_graph->create_initial_graph_multi_segment(_graph, initial_weight);
    }

    // Initialize aligner parameters and data structures
    _internal_penalties = InternalPenalties(penalties);
    _scope = std::make_unique<Scope>(n_scores);
    _beyond_scope = std::make_unique<BeyondScope>();
    constexpr int expected_nvertices = std::max(1024, 0); // TODO: Set the expected number of vertices
    _vertices_data = std::make_unique<VerticesData>(penalties, n_scores, expected_nvertices);
    _scratchpad = std::make_unique<ScratchPad>(-1024, 1024);
}

// Topological scope cache.
//
// Scope pruning for custom-end alignments previously ran a per-call BFS to
// find nodes that can reach end. That is O(graph size): a reverse BFS from a
// mid-graph end floods all its ancestors (~half the graph), even when start
// and end are adjacent anchors. Instead we cache a topological index per node
// (abPOA index_to_node_id pattern) and reuse it across many align_from calls.
void TheseusAlignerImpl::build_topo() {
    // Kahn's algorithm: assign each node a topological index. One linear pass
    // over the whole graph, amortised across many align_from calls.
    const size_t n = _graph.node_id_bound();
    _topo_pos.assign(n, -1);
    _topo_max_diag = 0;

    std::vector<int> indeg(n, 0);
    for (NodeId id : _graph.nodes()) {
        NodeView nv = _graph.node(id);
        for (auto out_id : nv.out_nodes)
            if (static_cast<size_t>(out_id) < n) ++indeg[out_id];
    }

    std::queue<NodeId> q;
    for (NodeId id : _graph.nodes()) {
        if (indeg[id] == 0) q.push(id);
        _topo_max_diag = std::max(_topo_max_diag, _graph.node_size(id));
    }

    int32_t idx = 0;
    while (!q.empty()) {
        NodeId cur = q.front(); q.pop();
        _topo_pos[cur] = idx++;
        NodeView nv = _graph.node(cur);
        for (auto out_id : nv.out_nodes) {
            if (static_cast<size_t>(out_id) < n && --indeg[out_id] == 0)
                q.push(out_id);
        }
    }
    // Any node left at -1 sits on a cycle; in_subgraph treats -1 as permissive.
    _topo_built_bound = static_cast<int>(n);
}

void TheseusAlignerImpl::maybe_build_topo() {
    // Rebuild when the graph has grown past the last build. Cheap amortisation:
    // most align_from calls add only a handful of nodes, so a build serves many
    // calls. Threshold avoids rebuilding for every single inserted base.
    const int bound = static_cast<int>(_graph.node_id_bound());
    if (_topo_built_bound < 0 || bound > _topo_built_bound + 32) {
        build_topo();
    }
}

void TheseusAlignerImpl::compute_subgraph(NodeId /*start*/, NodeId end) {
    maybe_build_topo();
    if (static_cast<size_t>(end) < _topo_pos.size() && _topo_pos[end] >= 0)
        _end_topo_pos = _topo_pos[end];
    else
        _end_topo_pos = 0x7fffffff;   // unknown end: do not prune
}

// Get the node/reversed node depending on the alignment configuration
NodeView TheseusAlignerImpl::get_node(NodeId id) {
  if (!_reversed_alignment) return _graph.node(id);
  else return _graph.node_rev(id);
}

bool TheseusAlignerImpl::has_out_nodes(NodeId id) {
  if (!_reversed_alignment) return !_graph.is_sink(id);
  else return !_graph.is_source(id);
}

void TheseusAlignerImpl::new_alignment(SequenceView seq,
                                       bool reverse_alignment,
                                       bool density_drop_active,
                                       bool lag_pruning_active) {
    _scope->new_alignment();
    _beyond_scope->new_alignment();
    _vertices_data->new_alignment();
    _seq = seq;
    _reversed_alignment = reverse_alignment;
    // Initialize scratchpad — only consider subgraph nodes if scoped.
    int max_diag = 0, v_n;
    if (_use_subgraph) {
        // max node size was accumulated during build_topo, avoiding an
        // O(graph size) scan here.
        max_diag = _topo_max_diag;
    } else {
        NodeIdRange nodes = _graph.nodes();
        for (const auto &node : nodes) {
            v_n = _graph.node_size(node);
            max_diag = std::max(max_diag, v_n);
        }
    }
    const int min_diag = -_seq.size();
    if (_scratchpad->max_diag() < max_diag ||
        _scratchpad->min_diag() > min_diag) {
        // TODO: Compute the max and min with a factor.
        _scratchpad = std::make_unique<ScratchPad>(min_diag, max_diag);
    }
    // Set data for first score
    _scope->new_score(_score);
    // Set initial alignment status
    _alignment.theseus_status = THESEUS_STATUS_OK;
    // Set heuristics
    _heuristics.new_alignment(_internal_penalties.gape(), _seq.size(), density_drop_active, lag_pruning_active);
    // TODO: Allow for different initial conditions. Now only global alignment.
    Cell init_condition;
    init_condition.offset    = 0;
    init_condition.vertex_id = _start_node;
    init_condition.diag = _start_offset;
    init_condition.prev_pos = -1;
    // Initial vertex data
    _beyond_scope->m_jumps_wf().push_back(init_condition);
    _vertices_data->activate_vertex(_start_node);
    _vertices_data->get_vertex_data(_start_node)._m_jumps_positions[0].push_back(0);
    // Alignment data
    _alignment.path.clear();
    _alignment.edit_op.clear();
}


// Process a given vertex with a given _score
void TheseusAlignerImpl::process_vertex(NodeId curr_node_id) {

  // Next
  int upper_bound = _graph.node_size(curr_node_id);
  next_I(upper_bound, curr_node_id);
  _scratchpad->reset();
  next_D(upper_bound, curr_node_id);
  _scratchpad->reset();
  next_M(upper_bound, curr_node_id);
  _scratchpad->reset();
  // Extend
  int v_pos = _vertices_data->get_id(curr_node_id);
  Scope::range cells_range = _scope->m_pos(_score)[v_pos];
  for (Cell::pos_t idx = cells_range.start; idx < cells_range.end; ++idx) {
    extend_diagonal(curr_node_id, _beyond_scope->m_wf()[idx], _beyond_scope->m_wf()[idx], idx, Cell::Matrix::M);
  }
}


void TheseusAlignerImpl::compute_new_wave() {
  // Update invalid segments
  _vertices_data->expand();
  _vertices_data->compact();
  // Process all active vertices
  int num_active_vertices = _vertices_data->num_active_vertices();
  NodeId curr_node_id;
  for (int l = 0; l < num_active_vertices; ++l) {
    curr_node_id = _vertices_data->get_vertex_id(l);
    process_vertex(curr_node_id);
  }
  _beyond_scope->update_positions();
}


Alignment TheseusAlignerImpl::align(
    std::string_view seq,
    int  weight,
    bool reverse_alignment,
    bool is_ends_free,
    int custom_start_node,
    int custom_start_offset,
    int custom_end_node,
    bool density_drop_active,
    bool lag_pruning_active
  )
{
  bool _custom_start = (custom_start_node >= 0);
  if (_custom_start) {
    _start_node = static_cast<NodeId>(custom_start_node);
    _start_offset = custom_start_offset;
    _end_node = (custom_end_node >= 0) ? static_cast<NodeId>(custom_end_node) : _sink_node;
  }
  else if (!reverse_alignment) {
    _start_node = _source_node;
    _start_offset = 0;
    _end_node = _sink_node;
  }
  // Align a sequence to a reference graph starting at a given position
  else {
    _start_node = _sink_node;
    _start_offset = 0;
    _end_node = _source_node;
  }
  _ends_free = is_ends_free;

  auto t0 = std::chrono::steady_clock::now();

  // Compute subgraph scope when a custom end node is set.
  if (_custom_start && custom_end_node >= 0) {
    _use_subgraph = true;
    compute_subgraph(_start_node, _end_node);
  } else {
    _use_subgraph = false;
  }

  auto t1 = std::chrono::steady_clock::now();

  // Initialize data for the new alignment
  new_alignment(Graph::SequenceView(seq, reverse_alignment), reverse_alignment, density_drop_active, lag_pruning_active);

  auto t2 = std::chrono::steady_clock::now();

  _score = 0;
  // _graph.print_code_graphviz();
  // Main alignment loop
  while (_alignment.theseus_status == THESEUS_STATUS_OK) {
    // Initial extend
    if (_score == 0) {
      NodeView start_node = get_node(_start_node);
      extend_diagonal(_start_node, _beyond_scope->m_jumps_wf()[0], _beyond_scope->m_jumps_wf()[0], 0, Cell::Matrix::MJumps);
    }
    // Next wave + extend
    compute_new_wave();
    // Evaluate global heuristics
    _alignment.theseus_status = (_alignment.theseus_status == THESEUS_STATUS_ALG_COMPLETED) ?
                                 _alignment.theseus_status :
                                  _heuristics.check_global_heuristics(_score);
    // Update score
    _score = _score + 1;
    // Clear the corresponding waves and metadata from the scope
    _scope->new_score(_score);
    _vertices_data->new_score(_score);
  }
  _score -= 1;

  auto t3 = std::chrono::steady_clock::now();

  // Backtrace
  _seq_ID += 1;
  if (_alignment.theseus_status == THESEUS_STATUS_ALG_COMPLETED) {
    backtrace();
    int start_column = _reversed_alignment ?
                       _graph.node_size(_start_pos.vertex_id) - (_start_pos.offset + _start_pos.diag) :
                       _start_pos.offset + _start_pos.diag;
    int start_row = _reversed_alignment ?
                    _seq.size() - _start_pos.offset :
                    0;
    if (_custom_start) {
      // Custom start: determine end-to-end by checking whether the path
      // actually reached the sink, since the alignment started at an
      // interior node (not the source).
      bool reached_end = false;
      if (!_alignment.path.empty()) {
          NodeId last_node = _reversed_alignment ? _alignment.path.front() : _alignment.path.back();
          reached_end = (last_node == _end_node);
      }
      _poa_graph->add_alignment_poa(_graph, _alignment, _seq, _seq_ID, start_column, start_row, weight, reached_end, _reversed_alignment, true, _start_offset);
    }
    else {
      // Original pericles path — unchanged
      _poa_graph->add_alignment_poa(_graph, _alignment, _seq, _seq_ID, start_column, start_row, weight, !_ends_free, _reversed_alignment);
    }
  }
  else {
    // Error messages
    if (_alignment.theseus_status == THESEUS_STATUS_MAX_STEPS_REACHED) {
      std::cerr << THESEUS_STATUS_MAX_STEPS_REACHED_MSG << std::endl;
    }
    else if (_alignment.theseus_status == THESEUS_STATUS_END_UNREACHABLE) {
      std::cerr << THESEUS_STATUS_END_UNREACHABLE_MSG << std::endl;
      // Choose starting position for backtrace (cell with maximum offset in M[s-_s_min]-M[s-_s_min - scope_size])
      // matrix on the last scope scores
      init_partial_backtrace();
      // Perform partial backtrace
      backtrace();
      int start_column = _reversed_alignment ?
                       _graph.node_size(_start_pos.vertex_id) - (_start_pos.offset + _start_pos.diag) :
                       _start_pos.offset + _start_pos.diag;
      int start_row = _reversed_alignment ?
                    _seq.size() - _start_pos.offset :
                    0;
      if (_custom_start) {
        _poa_graph->add_alignment_poa(_graph, _alignment, _seq, _seq_ID, start_column, start_row, weight, false, _reversed_alignment, true, _start_offset);
      }
      else {
        // Original pericles path — unchanged
        _poa_graph->add_alignment_poa(_graph, _alignment, _seq, _seq_ID, start_column, start_row, weight, !_ends_free, _reversed_alignment);
      }
    }
  }
  auto t4 = std::chrono::steady_clock::now();

  g_subgraph_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
  g_new_alignment_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
  g_dp_loop_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count();
  g_poa_update_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();
  g_total_score += _score;
  g_align_calls++;

  return _alignment;
}

  // Sparsify M data
  void TheseusAlignerImpl::sparsify_M_data(Cell::CellVector & dense_wf,
                                           int offset_increase,
                                           int shift_factor,
                                           Scope::range cells_range,
                                           int m,
                                           int upper_bound)
  {
    Cell::pos_t len = cells_range.end - cells_range.start, new_col;
    Cell new_cell;
    // Sparsify the active diagonals
    for (int l = 0; l < len; ++l)
    {
      new_cell = dense_wf[cells_range.start + l];
      new_cell.diag += shift_factor;
      new_cell.offset += offset_increase;
      new_cell.from_matrix = Cell::Matrix::M;
      new_cell.prev_pos = cells_range.start + l;
      new_col = new_cell.offset + new_cell.diag; // d = j - i -> j = d + i
      // Check validity
      if (new_cell.offset <= m && new_col <= upper_bound)
      { // If in bounds
        // Branchless push_back
        auto &cell = _scratchpad->access_alloc(new_cell.diag);
        // If better offset
        const bool cmp = cell.offset < new_cell.offset;
        cell = (cmp) ? new_cell : cell;
      }
    }
  }

  // Sparsify jumps data
  void TheseusAlignerImpl::sparsify_jumps_data(Cell::CellVector & dense_wf,
                                               std::vector<Cell::pos_t> & jumps_positions,
                                               int offset_increase,
                                               int shift_factor,
                                               int m,
                                               int upper_bound,
                                               Cell::Matrix from_matrix)
  {
    int len = jumps_positions.size(), new_col, pos;
    Cell new_cell;
    // Sparsify the active diagonals
    for (int l = 0; l < len; ++l)
    {
      pos = jumps_positions[l];
      new_cell = dense_wf[pos];
      new_cell.prev_pos = pos;
      new_cell.from_matrix = from_matrix;
      new_cell.diag += shift_factor;
      new_cell.offset += offset_increase;
      new_col = new_cell.offset + new_cell.diag; // d = j - i -> j = d + i
      // Check validity
      if (new_cell.offset <= m && new_col <= upper_bound)
      { // If in bounds
        // Branchless push_back
        auto &cell = _scratchpad->access_alloc(new_cell.diag);
        // If better offset
        const bool cmp = cell.offset < new_cell.offset;
        cell = (cmp) ? new_cell : cell;
      }
    }
  }

  // Sparsify indel
  void TheseusAlignerImpl::sparsify_indel_data(Cell::CellVector & dense_wf,
                                               int offset_increase,
                                               int shift_factor,
                                               Scope::range cells_range,
                                               int m,
                                               int upper_bound)
  {
    Cell::pos_t len = cells_range.end - cells_range.start, new_col;
    Cell new_cell;
    // Sparsify the active diagonals
    for (int l = 0; l < len; ++l)
    {
      // Vertex_id and previous matrix are the same as before
      new_cell = dense_wf[cells_range.start + l];
      new_cell.diag += shift_factor;
      new_cell.offset += offset_increase;
      new_col = new_cell.offset + new_cell.diag; // d = j - i -> j = d + i
      // Check validity
      if (new_cell.offset <= m && new_col <= upper_bound)
      { // If in bounds
        // Branchless push_back
        auto &cell = _scratchpad->access_alloc(new_cell.diag);
        // If better offset
        const bool cmp = cell.offset < new_cell.offset;
        cell = (cmp) ? new_cell : cell;
      }
    }
  }

  // Compute next I matrix
  void TheseusAlignerImpl::next_I(int upper_bound,
                                  NodeId curr_node_id)
  {
    // Sparsify data (put it in the scratch pad)
    int pos_prev_M = _score - (_internal_penalties.gapo() + _internal_penalties.gape()), pos_prev_I = _score - _internal_penalties.gape(), pos_prev_M_scope = _vertices_data->get_pos(pos_prev_M);
    int pos_prev_I_scope = _vertices_data->get_pos(pos_prev_I);
    // Come from an Insertion
    if (pos_prev_I >= 0) {
      if (_scope->i_pos(pos_prev_I).size() > _vertices_data->get_id(curr_node_id))
      {
        Scope::range cells_range = _scope->i_pos(pos_prev_I)[_vertices_data->get_id(curr_node_id)];
        sparsify_indel_data(_scope->i_wf(pos_prev_I), 0, 1, cells_range, _seq.size(), upper_bound); // Sparsify I data
      };
      sparsify_jumps_data(
        _beyond_scope->i_jumps_wf(),
        _vertices_data->get_vertex_data(curr_node_id)._i_jumps_positions[pos_prev_I_scope],
        0, 1, _seq.size(), upper_bound, Cell::Matrix::IJumps);
    }
    // Come from M
    if (pos_prev_M >= 0) {
      if (_scope->m_pos(pos_prev_M).size() > _vertices_data->get_id(curr_node_id)) {
        Scope::range cells_range = _scope->m_pos(pos_prev_M)[_vertices_data->get_id(curr_node_id)];
        sparsify_M_data(_beyond_scope->m_wf(), 0, 1, cells_range, _seq.size(), upper_bound); // Sparsify M data
      }
      sparsify_jumps_data(_beyond_scope->m_jumps_wf(),
          _vertices_data->get_vertex_data(curr_node_id)._m_jumps_positions[pos_prev_M_scope],
          0, 1, _seq.size(), upper_bound, Cell::Matrix::MJumps);
    }
    // Densify data (store it in the big wavefront)
    Scope::range new_range;
    new_range.start = _scope->i_wf(_score).size();
    for (auto diag : _scratchpad->active_diags()) {
      if (_vertices_data->valid_diagonal<Cell::Matrix::I>(curr_node_id, diag)
      && !_heuristics.check_local_heuristics((*_scratchpad)[diag].offset)) {
        _scope->i_wf(_score).push_back((*_scratchpad)[diag]);     // Store Cell
      }
    }
    new_range.end = _scope->i_wf(_score).size();
    _scope->i_pos(_score).push_back(new_range);
    // Check, store and invalidate new I jumps
    if (has_out_nodes(curr_node_id)) {
      NodeView curr_node = get_node(curr_node_id);
      check_and_store_jumps(curr_node, _scope->i_wf(_score), new_range);
    }
}


// Compute next D matrix
void TheseusAlignerImpl::next_D(int upper_bound,
                                NodeId curr_node_id)
{
  // Sparsify data (put it in the scratch pad)
  int pos_prev_M = _score - (_internal_penalties.gapo() + _internal_penalties.gape()),
      pos_prev_D = _score - _internal_penalties.gape(),
      pos_prev_M_scope = _vertices_data->get_pos(pos_prev_M);
  // Come from a Deletion
  if (pos_prev_D >= 0 && _scope->d_pos(pos_prev_D).size() > _vertices_data->get_id(curr_node_id))
  {
    Scope::range cells_range = _scope->d_pos(pos_prev_D)[_vertices_data->get_id(curr_node_id)];
    sparsify_indel_data(_scope->d_wf(pos_prev_D), 1, -1, cells_range, _seq.size(), upper_bound);
  }
  // Come from M
  if (pos_prev_M >= 0) {
    if (_scope->m_pos(pos_prev_M).size() > _vertices_data->get_id(curr_node_id))
    {
      Scope::range cells_range = _scope->m_pos(pos_prev_M)[_vertices_data->get_id(curr_node_id)];
      sparsify_M_data(_beyond_scope->m_wf(), 1, -1, cells_range, _seq.size(), upper_bound); // Sparsify M data
    }
    sparsify_jumps_data(_beyond_scope->m_jumps_wf(),
        _vertices_data->get_vertex_data(curr_node_id)._m_jumps_positions[pos_prev_M_scope],
        1, -1, _seq.size(), upper_bound, Cell::Matrix::MJumps);
  }
  // Densify data (store it in the big wavefront)
  Scope::range new_range;
  new_range.start = _scope->d_wf(_score).size();
  for (auto diag : _scratchpad->active_diags()) {
    if (_vertices_data->valid_diagonal<Cell::Matrix::D>(curr_node_id, diag)
    && !_heuristics.check_local_heuristics((*_scratchpad)[diag].offset)) {
      _scope->d_wf(_score).push_back((*_scratchpad)[diag]); // Store Cell
    }
  }
  new_range.end = _scope->d_wf(_score).size();
  _scope->d_pos(_score).push_back(new_range);
}


// Compute next M matrix
void TheseusAlignerImpl::next_M(int upper_bound,
                                NodeId curr_node_id) {
  // Sparsify data (put it in the scratch pad)
  int pos_prev_M = _score - _internal_penalties.mism(),
      pos_prev_D = _score,
      pos_prev_I = _score,
      pos_prev_M_scope = _vertices_data->get_pos(pos_prev_M);
  // Come from a Deletion
  if (_scope->d_pos(pos_prev_D).size() > _vertices_data->get_id(curr_node_id))  {
    Scope::range cells_range = _scope->d_pos(pos_prev_D)[_vertices_data->get_id(curr_node_id)];
    sparsify_indel_data(_scope->d_wf(pos_prev_D), 0, 0, cells_range, _seq.size(), upper_bound);
  }
  // Come from an Insertion
  if (_scope->i_pos(pos_prev_I).size() > _vertices_data->get_id(curr_node_id))  {
    Scope::range cells_range = _scope->i_pos(pos_prev_I)[_vertices_data->get_id(curr_node_id)];
    sparsify_indel_data(_scope->i_wf(pos_prev_I), 0, 0, cells_range, _seq.size(), upper_bound);
  }
  // Come from M
  if (pos_prev_M >= 0) {
    if (_scope->m_pos(pos_prev_M).size() > _vertices_data->get_id(curr_node_id))  {
      Scope::range cells_range = _scope->m_pos(pos_prev_M)[_vertices_data->get_id(curr_node_id)];
      sparsify_M_data(_beyond_scope->m_wf(), 1, 0, cells_range,  _seq.size(), upper_bound);
    }
    sparsify_jumps_data(_beyond_scope->m_jumps_wf(),
        _vertices_data->get_vertex_data(curr_node_id)._m_jumps_positions[pos_prev_M_scope],
        1, 0, _seq.size(), upper_bound, Cell::Matrix::MJumps);
  }
  // Densify data (store it in the big wavefront)
  Scope::range new_range;
  new_range.start = _beyond_scope->m_wf().size();
  for (auto diag : _scratchpad->active_diags()) {
    if (_vertices_data->valid_diagonal<Cell::Matrix::M>(curr_node_id, diag)
    && !_heuristics.check_local_heuristics((*_scratchpad)[diag].offset)) {
      _beyond_scope->m_wf().push_back((*_scratchpad)[diag]);     // Store Cell
    }
  }
  new_range.end = _beyond_scope->m_wf().size();
  _scope->m_pos(_score).push_back(new_range);
}


// Store the jump in neighbours
void TheseusAlignerImpl::store_M_jump(NodeView curr_node,
                                      Cell &prev_cell,
                                      Cell::pos_t prev_pos,
                                      Cell::Matrix from_matrix) {
  // Invalidate the jumping diagonal
  // if (_seq_ID >= 2653 && prev_cell.vertex_id == 351) {
  //    std::cout << "Invalidating diag " << prev_cell.diag << " in vertex " << prev_cell.vertex_id << std::endl;
  //    int vertex_id = _vertices_data->get_id(prev_cell.vertex_id);
  // }
  _vertices_data->invalidate_m_jump(_vertices_data->get_id(prev_cell.vertex_id), prev_cell.diag);
  int pos_score = _vertices_data->get_pos(_score);
  int new_diag  = -prev_cell.offset;
  Cell new_cell = prev_cell;
  new_cell.from_matrix = from_matrix;
  new_cell.prev_pos = prev_pos;
  // For each neighbour, store the jump and metadata.
  // Skip neighbours outside the subgraph when scoped.
  for (auto out_node_id : curr_node.out_nodes) {
    if (!in_subgraph(out_node_id)) continue;
    new_cell.vertex_id = out_node_id;
    new_cell.diag = new_diag;
    NodeView out_node = get_node(out_node_id);
    _vertices_data->activate_vertex(new_cell.vertex_id);
    // Store jump and metadata
    bool valid_diag = _vertices_data->valid_diagonal<Cell::Matrix::M>(new_cell.vertex_id, new_cell.diag);
    // Extend only if it has not yet been visited
    if (valid_diag) {
      int pos_new_cell = _beyond_scope->m_jumps_wf().size();
      _beyond_scope->m_jumps_wf().push_back(new_cell);
      _vertices_data->get_vertex_data(new_cell.vertex_id)._m_jumps_positions[pos_score].push_back(pos_new_cell);
      extend_diagonal(out_node_id, _beyond_scope->m_jumps_wf()[pos_new_cell], _beyond_scope->m_jumps_wf()[pos_new_cell], pos_new_cell, Cell::Matrix::MJumps);
    }
  }
}


// Store the jump in neighbours
void TheseusAlignerImpl::store_I_jump(
    NodeView curr_node,
    Cell& prev_cell,
    Cell::pos_t prev_pos,
    Cell::Matrix from_matrix)
{
  // Invalidate the jumping diagonal
  _vertices_data->invalidate_i_jump(_vertices_data->get_id(prev_cell.vertex_id), prev_cell.diag);
  int pos_score = _vertices_data->get_pos(_score);
  int new_diag = -prev_cell.offset;
  Cell new_cell = prev_cell;
  new_cell.from_matrix = from_matrix;
  new_cell.prev_pos = prev_pos;
  // For each neighbour, store the jump and metadata.
  // Skip neighbours outside the subgraph when scoped.
  for (auto out_node_id : curr_node.out_nodes) {
    if (!in_subgraph(out_node_id)) continue;
    new_cell.vertex_id = out_node_id;
    new_cell.diag = new_diag;
    _vertices_data->activate_vertex(new_cell.vertex_id);
    bool valid_diag = _vertices_data->valid_diagonal<Cell::Matrix::I>(new_cell.vertex_id, new_cell.diag);
    // Extend only if it has not yet been visited
    if (valid_diag) {
      int pos_new_cell = _beyond_scope->i_jumps_wf().size();
      _beyond_scope->i_jumps_wf().push_back(new_cell);
      _vertices_data->get_vertex_data(new_cell.vertex_id)._i_jumps_positions[pos_score].push_back(pos_new_cell);
      // If the destination vertex is empty, jump again
      if (curr_node.sequence.empty()) {
        store_I_jump(curr_node, _beyond_scope->i_jumps_wf()[pos_new_cell], prev_pos, Cell::Matrix::IJumps);
      }
    }
  }
}


// Check and store I jumps (that is, those diagonals that have reached the last column of a vertex)
void TheseusAlignerImpl::check_and_store_jumps(
    NodeView curr_node,
    Cell::CellVector &curr_wavefront,
    Scope::range cell_range)
{
  Cell::pos_t len = cell_range.end - cell_range.start, n = curr_node.sequence.size(), prev_pos;
  Cell::idx2d_t diag, offset, curr_j;
  Cell::Matrix from_matrix;
  // Check all diagonals in the current wavefront
  for (int l = 0; l < len; ++l) {
    diag = curr_wavefront[cell_range.start + l].diag;
    offset = curr_wavefront[cell_range.start + l].offset;
    curr_j = diag + offset;
    if (curr_j == n && offset <= (int)_seq.size()) {
      from_matrix = curr_wavefront[cell_range.start + l].from_matrix;
      prev_pos = curr_wavefront[cell_range.start + l].prev_pos;
      store_M_jump(curr_node, curr_wavefront[cell_range.start + l], prev_pos, from_matrix);
      store_I_jump(curr_node, curr_wavefront[cell_range.start + l], prev_pos, from_matrix);
    }
  }
}


// Compute the Longest Common Prefix between two given sequences
void TheseusAlignerImpl::LCP(
    NodeView &curr_node,
    int &offset,
    int &j)
{
  // Find LCP
  int len_seq_1 = _seq.size();
  int len_seq_2 = curr_node.sequence.size();
  // std::cout << query[offset] << " " << curr_node.sequence[j] << std::endl;
  while (offset < len_seq_1 && j < len_seq_2 && _seq[offset] == curr_node.sequence[j]) {
    offset = offset + 1;   // Update the f.r. of this diagonal
    j = j + 1;
  }
}


// TODO: Implement different end conditions as Global, Semi-Global...
void TheseusAlignerImpl::check_end_condition(
    Cell curr_data)
{
  // Ends free
  if (_ends_free) {
    if (curr_data.offset == (int)_seq.size() || curr_data.vertex_id == _end_node) {
    _alignment.theseus_status = THESEUS_STATUS_ALG_COMPLETED;
    _start_pos = curr_data;
    }
  }
  // End to end
  else {
    if (curr_data.offset == (int)_seq.size() && curr_data.vertex_id == _end_node) {
      _alignment.theseus_status = THESEUS_STATUS_ALG_COMPLETED;
      _start_pos = curr_data;
    }
  }
}


// Extend a particular diagonal
void TheseusAlignerImpl::extend_diagonal(
    NodeId curr_node_id,
    Cell &curr_cell,
    Cell &prev_cell,
    Cell::pos_t prev_pos,
    Cell::Matrix from_matrix)
{
  // Longest Common prefix
  NodeView curr_node = get_node(curr_node_id);
  int j = curr_cell.diag + curr_cell.offset;
  LCP(curr_node, curr_cell.offset, j);
  // End condition
  check_end_condition(curr_cell);
  // _graph.print_graph();
  if (j == (int)curr_node.sequence.size() && curr_cell.offset <= (int)_seq.size() && has_out_nodes(curr_node_id)) {
    store_M_jump(curr_node, prev_cell, prev_pos, from_matrix); // Store the jump in neighbours
  }
}


// Initialize the partial backtrace, by finding the starting cell for backtrace.
// We find this cell by checking all the active cells in scores (s-_s_min, ..., s-_s_min-_scope_size)
void TheseusAlignerImpl::init_partial_backtrace() {
  Cell best_cell;
  best_cell.offset = -1;
  // Iterate through the valid score range
  int start_score = _score - _heuristics.s_min() - _scope->size();
  start_score = std::max(0, start_score);
  // Iterate in the M structure for score s TODO: I_jumps?
  for (int s = start_score; s <= _score - _heuristics.s_min(); ++s) {
    // Check M wavefront
    int start_pos_M = (s == 0) ? 0 : _beyond_scope->m_wf_pos(s-1);
    int end_pos_M   = _beyond_scope->m_wf_pos(s) - 1;
    for (int pos = start_pos_M; pos < end_pos_M; ++pos) {
      Cell curr_cell = _beyond_scope->m_wf()[pos];
      if (curr_cell.offset > best_cell.offset) {
        best_cell = curr_cell;
      }
    }
    // Check M_jumps wavefront
    int start_pos_M_jumps = (s == 0) ? 0 : _beyond_scope->m_jumps_wf_pos(s-1);
    int end_pos_M_jumps   = _beyond_scope->m_jumps_wf_pos(s) - 1;
    for (int pos = start_pos_M_jumps; pos < end_pos_M_jumps; ++pos) {
      Cell curr_cell = _beyond_scope->m_jumps_wf()[pos];
      if (curr_cell.offset > best_cell.offset) {
        best_cell = curr_cell;
      }
    }
  }
  _start_pos = best_cell;
}


// Add matches to our backtracking vector
void TheseusAlignerImpl::add_matches(
    int start_matches,
    int end_matches)
{
  int size = end_matches - start_matches;
  for (int k = 0; k < size; ++k) {
    _alignment.edit_op.push_back('M');
  }
}


// Add a mismatch to our backtracking vector
void TheseusAlignerImpl::add_mismatch()
{
  _alignment.edit_op.push_back('X');
}


// Add an insertion to our backtracking vector
void TheseusAlignerImpl::add_insertion()
{
  _alignment.edit_op.push_back('I');
}


// Add a deletion to our backtracking vector
void TheseusAlignerImpl::add_deletion()
{
  _alignment.edit_op.push_back('D');
}


void TheseusAlignerImpl::one_backtrace_step(
    Cell &curr_cell)
{
  Cell prev_cell;
  if (curr_cell.from_matrix == Cell::Matrix::M) prev_cell = _beyond_scope->m_wf()[curr_cell.prev_pos];
  else if (curr_cell.from_matrix == Cell::Matrix::MJumps) prev_cell = _beyond_scope->m_jumps_wf()[curr_cell.prev_pos];
  else prev_cell = _beyond_scope->i_jumps_wf()[curr_cell.prev_pos];
  int num_indels;
  bool is_jump = ((curr_cell.vertex_id != prev_cell.vertex_id));
  // We are inside the same vertex
  if (!is_jump) { // Still in the same vertex
    // Mismatch
    if (curr_cell.diag == prev_cell.diag) {
      if (curr_cell.offset > prev_cell.offset) {    // Consider 0 length vertices
        add_matches(prev_cell.offset + 1, curr_cell.offset);
        add_mismatch();
      }
    }
    else {
      // Deletion
      if (curr_cell.diag < prev_cell.diag) {
        num_indels = prev_cell.diag - curr_cell.diag;
        add_matches(prev_cell.offset + num_indels, curr_cell.offset);
        for (int l = 0; l < num_indels; ++l) add_deletion();
      }
      // Insertion
      else {
        num_indels = curr_cell.diag - prev_cell.diag;
        add_matches(prev_cell.offset, curr_cell.offset);
        for (int l = 0; l < num_indels; ++l) add_insertion();
      }
    }
  }
  // Jump to another vertex
  else {
    add_matches(prev_cell.offset, curr_cell.offset);                          // Add the necessary matches
    _alignment.path.push_back(prev_cell.vertex_id);                           // Add the new vertex to the path
    int col_in_prev_v = prev_cell.diag + prev_cell.offset;
    int num_insertions = _graph.node_size(prev_cell.vertex_id) - col_in_prev_v;
    for (int l = 0; l < num_insertions; ++l) add_insertion();                 // Add the necessary insertions
  }
  // Update current cell
  curr_cell = prev_cell;
}


// Main function of the backtracking process
void TheseusAlignerImpl::backtrace()
{
  Cell curr_pos = _start_pos;
  _alignment.start_offset = _start_offset;
  _alignment.end_offset = curr_pos.diag + curr_pos.offset; // Vertex offset = j
  _alignment.path.push_back(curr_pos.vertex_id);
  // Main backtrace loop
  while (curr_pos.prev_pos != -1)
  {
    one_backtrace_step(curr_pos);
  }
  // Add the matches until the beginning of the sequence
  add_matches(0, curr_pos.offset);
  // Reverse path and edit operations
  if (!_reversed_alignment) {
    std::reverse(_alignment.edit_op.begin(), _alignment.edit_op.end());
    std::reverse(_alignment.path.begin(), _alignment.path.end());
  }
}


// Output functions
/**
 * @brief Visualize the graph in Graphviz format.
 *
 * @param G
 */
void TheseusAlignerImpl::print_code_graphviz_internal(std::ostream &out_stream)
{
    // TODO: What if we wanted to print reversed?
    _reversed_alignment = false;
    out_stream << "digraph G {" << std::endl;
    // Print nodes
    for (NodeId id : _graph.nodes())
    {
      NodeView node = get_node(id);
      out_stream << id << " [label=\"";
      for (size_t j = 0; j < node.sequence.size(); ++j)
      {
          out_stream << node.sequence[j];
      }
      out_stream << "\"]" << std::endl;
    }
    // Print edges
    for (NodeId id : _graph.nodes())
    {
      auto node = get_node(id);
      for (NodeId out_id : node.out_nodes)
      {
        out_stream << id << "->" << out_id << std::endl;
      }
    }
    out_stream << "}" << std::endl;
}
void TheseusAlignerImpl::print_code_graphviz(std::ostream &out_stream) {
  print_code_graphviz_internal(out_stream);
}

// Print as GFA
void TheseusAlignerImpl::print_as_gfa_internal(
  std::ostream &gfa_output)
{
    // Print all nodes as Segments
    for (NodeId id : _graph.nodes())
    {
      NodeView node = get_node(id);
      gfa_output << "S\t" << id << "\t" << node.sequence << "\n";
    }
    // Print all edges as Links
    for (NodeId id : _graph.nodes())
    {
      NodeView node = get_node(id);
      // Go through all incoming vertices (with this you cover all possible edges,
      // since the graph is directed)
      for (const auto in_node : node.in_nodes)
      {
        gfa_output << "L\t" << in_node << "\t+\t"
          << id << "\t+\t" << "0M\n";
      }
    }
}
// TODO: Using node names?
void TheseusAlignerImpl::print_as_gfa(std::ostream &gfa_output)
{
  print_as_gfa_internal(gfa_output);
}

std::vector<uint8_t> TheseusAlignerImpl::get_msa_matrix(int num_sequences,
                                                        int &n_rows, int &n_cols) {
  return _poa_graph->poa_to_msa_matrix(
      num_sequences >= 0 ? num_sequences : _seq_ID, n_rows, n_cols);
}

// Print as msa (can only call from TheseusMSA)
void TheseusAlignerImpl::print_as_msa(std::ostream &out_stream, int num_sequences,
                                      const std::vector<std::string> *seq_names) {
  _poa_graph->poa_to_fasta(num_sequences >= 0 ? num_sequences : _seq_ID, out_stream, seq_names);
}


// Find and return the consensus sequence (can only call from TheseusMSA)
std::string TheseusAlignerImpl::heaviest_bundle_consensus() {
  return _poa_graph->poa_to_consensus();
}

// Find and return the consensus sequence using weighted majority voting
void TheseusAlignerImpl::majority_voting_consensus(std::vector<int> &consensus_weights,
                                                   std::string &consensus_sequence,
                                                   std::string &consensus_sequence_gapped) {
  _poa_graph->poa_to_consensus_weighted_majority_voting(_seq_ID, consensus_weights, consensus_sequence, consensus_sequence_gapped);
}

// Print as GAF
void TheseusAlignerImpl::print_as_gaf(
    theseus::Alignment &alignment,
    std::ostream &out_stream,
    std::string seq_name,
    std::unordered_map<NodeId, std::string> &node_names) {

  // Field 1: Query name
  out_stream << seq_name;

  // Field 2: Query length
  out_stream << "\t" << _seq.size();

  // Field 3: Query start
  out_stream << "\t" << 0;

  // Field 4: Query end
  out_stream << "\t" << _seq.size();

  // Field 5: Strand
  out_stream << "\t" << "+"; // TODO: Support reverse strand

  // Field 6: Alignment path
  out_stream << "\t";
  for (size_t l = 0; l < alignment.path.size(); ++l) {
    out_stream << ">" << node_names.at(alignment.path[l]); // TODO: Support orientation
  }

  // Field 7: Target length
  int target_length = 0;
  for (size_t l = 0; l < alignment.path.size(); ++l) {
    int node_length = node_names.at(alignment.path[l]).size();
    target_length += node_length;
  }
  out_stream << "\t" << target_length;

  // Field 8: Target start
  out_stream << "\t" << _start_offset;

  // Field 9: Target end
  out_stream << "\t" << alignment.end_offset;

  // Field 10: Number of matching bases
  int num_matches = 0;
  for (size_t l = 0; l < alignment.edit_op.size(); ++l) {
    if (alignment.edit_op[l] == 'M') {
      num_matches += 1;
    }
  }
  out_stream << "\t" << num_matches;

  // Field 11: Alignment block length
  out_stream << "\t" << alignment.edit_op.size();

  // Field 12: Mapping quality
  out_stream << "\t" << 255; // TODO: Compute mapping quality

  // Optional fields
  out_stream << "\t" << "cg:Z:"; // CIGAR string
  std::string cigar = "";
  int count = 1;
  for (size_t l = 1; l < alignment.edit_op.size(); ++l) {
    if (alignment.edit_op[l] == alignment.edit_op[l - 1]) {
      count += 1;
    }
    else {
      cigar += std::to_string(count) + alignment.edit_op[l - 1];
      count = 1;
    }
  }
  if (alignment.edit_op.size() > 0) {
    cigar += std::to_string(count) + alignment.edit_op.back();
  }
  out_stream << cigar << "\n";
}

} // namespace theseus