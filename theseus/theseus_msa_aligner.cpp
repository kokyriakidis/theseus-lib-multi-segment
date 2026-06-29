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


#include "theseus/theseus_msa_aligner.h"

#include "theseus_aligner_impl.h"
#include <atomic>

namespace theseus {

using NodeId = Graph::NodeId;

TheseusMSA::TheseusMSA(const Penalties &penalties,
                       const Heuristics &heuristics,
                       std::string_view seq,
                       int initial_weight,
                       bool is_ends_free) {

    // Error out if the sequence is ends_free. Ask for a backbone (end-to-end)
    if (is_ends_free) {
        throw std::invalid_argument("Sequence must be end-to-end for MSA initialization");
    }
    // Create the initial graph
    theseus::Graph G;
    // Add nodes
    NodeId source_node_id = G.add_node("");
    NodeId central_node_id = G.add_node(seq);
    NodeId sink_node_id = G.add_node("");
    // Add edges
    G.add_edge(source_node_id, central_node_id);
    G.add_edge(central_node_id, sink_node_id);
    // Construct the aligner implementation
    msa_aligner_impl_ = std::make_unique<TheseusAlignerImpl>(penalties, heuristics, std::move(G), initial_weight);
}

TheseusMSA::TheseusMSA(const Penalties &penalties,
                       const Heuristics &heuristics,
                       const std::vector<std::string_view> &segments,
                       std::vector<Graph::NodeId> &node_ids,
                       int initial_weight) {

    theseus::Graph G;
    // Source node (empty)
    NodeId source_node_id = G.add_node("");
    // One node per segment
    node_ids.clear();
    node_ids.reserve(segments.size());
    NodeId prev_node_id = source_node_id;
    for (const auto &seg : segments) {
        NodeId seg_node_id = G.add_node(seg);
        G.add_edge(prev_node_id, seg_node_id);
        node_ids.push_back(seg_node_id);
        prev_node_id = seg_node_id;
    }
    // Sink node (empty)
    NodeId sink_node_id = G.add_node("");
    G.add_edge(prev_node_id, sink_node_id);
    // Construct the aligner implementation with the correct sink node
    msa_aligner_impl_ = std::make_unique<TheseusAlignerImpl>(
        penalties, heuristics, std::move(G), initial_weight, sink_node_id);
}

/**
 * @brief Destroy the Theseus Aligner:: Theseus Aligner object
 *
 */
TheseusMSA::~TheseusMSA() {}

/**
 * @brief Main alignment function for the Theseus aligner.
 *
 * @param seq
 * @param weight              Sequence weight
 * @param reverse_alignment   Whether to align in reverse
 * @param is_ends_free        Whether to allow a free end
 * @param density_drop_active Whether to use the density-drop heuristic
 * @param lag_pruning_active  Whether to use the lag pruning heuristic
 * @return Alignment
 */
Alignment TheseusMSA::align(
    std::string_view seq,
    int weight,
    bool reverse_alignment,
    bool is_ends_free,
    bool density_drop_active,
    bool lag_pruning_active) {

    // Error out if you are adding a backbone sequence after partial sequences
    if (!still_end_to_end && !is_ends_free) {
        throw std::invalid_argument("Cannot add a backbone sequence after partial sequences");
    }
    // Update the still_end_to_end value
    still_end_to_end = !is_ends_free;
    return msa_aligner_impl_->align(seq, weight, reverse_alignment, is_ends_free,
                                    -1, 0, -1,
                                    density_drop_active, lag_pruning_active);
}

Alignment TheseusMSA::align_from(
    std::string_view seq,
    Graph::NodeId start_node,
    int weight,
    bool is_ends_free,
    int start_offset,
    int end_node,
    int seq_id) {

    // Partial alignments are always ends-free by nature
    still_end_to_end = false;
    // Override the next auto-assigned seq_id so that multiple segments
    // of the same read share one MSA row.
    if (seq_id >= 0) {
        msa_aligner_impl_->set_next_seq_id(seq_id);
    }
    return msa_aligner_impl_->align(seq, weight, false, is_ends_free,
                                    static_cast<int>(start_node), start_offset,
                                    end_node);
}

/**
 * @brief Print the current POA graph in MSA format.
 *
 */
void TheseusMSA::print_as_gfa(std::ostream &out_stream) {
    msa_aligner_impl_->print_as_gfa(out_stream);
}

/**
 * @brief Print the current POA graph in MSA format.
 *
 */
void TheseusMSA::print_as_msa(std::ostream &out_stream, int num_sequences,
                              const std::vector<std::string> *seq_names) {
    msa_aligner_impl_->print_as_msa(out_stream, num_sequences, seq_names);
}

MSAMatrix TheseusMSA::get_msa_matrix(int num_sequences,
                                     const std::vector<uint64_t> *seq_ids) {
    MSAMatrix m;
    m.data = msa_aligner_impl_->get_msa_matrix(num_sequences, m.n_rows, m.n_cols);
    if (seq_ids) {
        m.seq_ids = *seq_ids;
    }
    return m;
}

/**
 * @brief Return consensus sequence.
 *
 */
std::string TheseusMSA::heaviest_bundle_consensus() {
    return msa_aligner_impl_->heaviest_bundle_consensus();
}

/**
 * @brief Print the weighted majority voting consensus.
 *
 */
void TheseusMSA::majority_voting_consensus(std::vector<int> &consensus_weights,
                                           std::string &consensus_sequence,
                                           std::string &consensus_sequence_gapped) {
    msa_aligner_impl_->majority_voting_consensus(consensus_weights, consensus_sequence, consensus_sequence_gapped);
}

/**
 * @brief Print in graphviz format.
 *
 */
void TheseusMSA::print_as_dot(std::ostream &out_stream) {
    msa_aligner_impl_->print_code_graphviz(out_stream);
}

} // namespace theseus

// Timing counter access — references the atomics in theseus_aligner_impl.cpp.
extern std::atomic<int64_t> g_subgraph_ns;
extern std::atomic<int64_t> g_new_alignment_ns;
extern std::atomic<int64_t> g_dp_loop_ns;
extern std::atomic<int64_t> g_backtrace_ns;
extern std::atomic<int64_t> g_poa_update_ns;
extern std::atomic<int64_t> g_align_calls;
extern std::atomic<int64_t> g_total_score;

theseus::TheseusMSA::TimingCounters theseus::TheseusMSA::get_timing_counters() {
    TimingCounters tc;
    tc.subgraph_ns = g_subgraph_ns.load();
    tc.new_alignment_ns = g_new_alignment_ns.load();
    tc.dp_loop_ns = g_dp_loop_ns.load();
    tc.poa_update_ns = g_poa_update_ns.load();
    tc.align_calls = g_align_calls.load();
    tc.total_score = g_total_score.load();
    return tc;
}

void theseus::TheseusMSA::reset_timing_counters() {
    g_subgraph_ns = 0;
    g_new_alignment_ns = 0;
    g_dp_loop_ns = 0;
    g_backtrace_ns = 0;
    g_poa_update_ns = 0;
    g_align_calls = 0;
    g_total_score = 0;
}