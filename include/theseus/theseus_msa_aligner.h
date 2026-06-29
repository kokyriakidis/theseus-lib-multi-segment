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


#pragma once

#include <memory>
#include <istream>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "theseus/penalties.h"
#include "theseus/alignment.h"
#include "theseus/heuristics.h"
#include "theseus/graph.h"

/**
 * Multiple Sequence Aligner (MSA) based on POA graphs. Internally uses the
 * TheseusAlignerImpl class.
 *
 */

namespace theseus
{

    /** Row-major MSA matrix with sequence IDs.
     *  Access element at (row, col) via data[row * n_cols + col].
     *  Values: 'A','C','G','T' for bases, '-' for gaps.
     *  seq_ids[i] is the caller-provided numeric ID for row i. */
    struct MSAMatrix {
        std::vector<uint8_t> data;       // flat row-major storage
        std::vector<uint64_t> seq_ids;   // one numeric ID per row
        int n_rows = 0;                  // number of sequences (including backbone)
        int n_cols = 0;                  // number of MSA columns

        uint8_t  operator()(int row, int col) const { return data[row * n_cols + col]; }
        uint8_t& operator()(int row, int col)       { return data[row * n_cols + col]; }
    };

    class TheseusAlignerImpl; // Forward declaration of the implementation class.

    class TheseusMSA
    {
    public:
        /**
         * Initial constructor
         *
         * @param penalties    User defined alignment penalties
         * @param heuristics   User defined heuristics
         * @param seq          Sequence to initialize the graph
         * @param initial_weight    Initial weight for the first added sequence
         * @param is_ends_free Whether to use ends-free alignment
         */
        TheseusMSA(
            const Penalties &penalties,
            const Heuristics &heuristics,
            std::string_view seq,
            int initial_weight,
            bool is_ends_free);

        /**
         * Multi-segment constructor. Creates a graph with one node per segment,
         * chained as: source -> seg0 -> seg1 -> ... -> segN -> sink.
         *
         * This is used when the backbone is split at anchor boundaries so that
         * each anchor position corresponds to a node boundary in the graph.
         *
         * @param penalties      User defined alignment penalties
         * @param heuristics     User defined heuristics
         * @param segments       Vector of sequences, one per backbone segment
         * @param node_ids       Output: filled with the compact graph NodeId for
         *                       each segment (node_ids[i] is the NodeId of segments[i])
         * @param initial_weight Initial weight for the backbone sequence
         */
        TheseusMSA(
            const Penalties &penalties,
            const Heuristics &heuristics,
            const std::vector<std::string_view> &segments,
            std::vector<Graph::NodeId> &node_ids,
            int initial_weight = 1);

        /**
         * Class destructor
         *
         */
        ~TheseusMSA();

        /**
         * Add a new sequence to the POA graph, representing the MSA so far.
         *
         * @param seq Sequence to add to the MSA
         */
        Alignment align(std::string_view seq,
                        int  weight = 1,
                        bool reverse_alignment = false,
                        bool is_ends_free = false,
                        bool density_drop_active = false,
                        bool lag_pruning_active = false);

        /**
         * Add a new sequence to the POA graph starting at a specific node.
         * The alignment is scoped to the subgraph between start_node and
         * end_node. Only nodes reachable from start_node that can reach
         * end_node are considered during alignment.
         *
         * @param seq           Sequence to add to the MSA
         * @param start_node    Node ID to start alignment from
         * @param weight        Weight of the sequence (default 1)
         * @param is_ends_free  Whether to allow a free end (default true)
         * @param start_offset  Offset within the starting node (default 0)
         * @param end_node      Node ID to end alignment at (-1 = sink, default)
         * @param seq_id        Sequence ID for MSA row assignment (-1 = auto-increment, default).
         *                      When >= 0, the internally auto-assigned ID is remapped to this
         *                      value in all POA vertices, allowing multiple segments to share
         *                      one MSA row.
         */
        Alignment align_from(std::string_view seq,
                             Graph::NodeId start_node,
                             int weight = 1,
                             bool is_ends_free = true,
                             int start_offset = 0,
                             int end_node = -1,
                             int seq_id = -1);

        /**
         * @brief Print the current POA graph as a GFA file.
         *
         */
        void print_as_gfa(std::ostream &out_stream);


        /**
         * @brief Print the current POA graph in MSA format.
         *
         * @param out_stream    Output stream
         * @param num_sequences Override the number of sequences (rows) in the MSA.
         *                      -1 (default) uses the internal auto-incremented count.
         *                      Use this when seq_id remapping has reduced the
         *                      number of distinct sequences.
         */
        void print_as_msa(std::ostream &out_stream, int num_sequences = -1,
                          const std::vector<std::string> *seq_names = nullptr);


        /**
         * @brief Compute and return the MSA matrix.
         *
         * @param num_sequences  Number of sequences (rows). -1 = auto.
         * @param seq_names      Sequence names for each row. If provided,
         *                       stored in the returned MSAMatrix.
         * @return MSAMatrix     Flat row-major matrix with base/gap values and names.
         */
        MSAMatrix get_msa_matrix(int num_sequences = -1,
                                 const std::vector<uint64_t> *seq_ids = nullptr);

        /**
         * @brief Return consensus sequence.
         *
         */
        std::string heaviest_bundle_consensus();

        /**
         * @brief Compute the weighted majority voting consensus sequence.
         *
         */
        void majority_voting_consensus(std::vector<int> &consensus_weights,
                                       std::string &consensus_sequence,
                                       std::string &consensus_sequence_gapped);


        /**
         * @brief Print in graphviz format.
         *
         */
        void print_as_dot(std::ostream &out_stream);

        /**
         * @brief Cumulative timing breakdown from all align/align_from calls.
         * Values are in nanoseconds.
         */
        struct TimingCounters {
            int64_t subgraph_ns = 0;
            int64_t new_alignment_ns = 0;
            int64_t dp_loop_ns = 0;
            int64_t poa_update_ns = 0;
            int64_t align_calls = 0;
            int64_t total_score = 0;
        };
        static TimingCounters get_timing_counters();
        static void reset_timing_counters();


    private:
        std::unique_ptr<TheseusAlignerImpl> msa_aligner_impl_;
        bool still_end_to_end = true;
    };

} // namespace theseus
