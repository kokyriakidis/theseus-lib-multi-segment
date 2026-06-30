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


#include "../doctest.h"

#include <vector>
#include <string>
#include <iostream>
#include "../../include/theseus/alignment.h"
#include "../../include/theseus/penalties.h"
#include "../../include/theseus/theseus_msa_aligner.h"
#include "../../include/theseus/heuristics.h"
#include "../../include/theseus/graph.h"

using NodeId = theseus::Graph::NodeId;


TEST_CASE("Check MSA aligner") {
    SUBCASE("Correct MSA with a matching sequence") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "ACCCGTAAAAGGG";
        std::vector<char> expected_cigar = {'M','M','M','M','M','M','M','M','M','M','M','M','M'};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1,false); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, 1, false);

        // Check if alignment was successful
        CHECK(alignment.compute_affine_gap_score(penalties) == 0); // Check score
        CHECK(alignment.edit_op == expected_cigar);                  // Check CIGAR
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));     // Check path
    }

    SUBCASE("Correct MSA with a sequence with a mismatch") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "ACCCGTCAAAGGG";
        std::vector<char> expected_cigar = {'M','M','M','M','M','M','X','M','M','M','M','M','M'};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1,false); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, 1, false);

        // Check if alignment was successful
        CHECK(alignment.compute_affine_gap_score(penalties) == 2); // Check score
        CHECK(alignment.edit_op == expected_cigar);                  // Check CIGAR
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));        // Check path
    }

    SUBCASE("Correct MSA with an deletion at the end") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "ACCCGTAAAAGGGAAA";
        std::vector<char> expected_cigar = {'M','M','M','M','M','M','M','M','M','M','M','M','M', 'D','D','D'};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1,false); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, false);

        // Check if alignment was successful
        CHECK(alignment.compute_affine_gap_score(penalties) == 6);    // Check score
        CHECK(alignment.edit_op == expected_cigar);                      // Check CIGAR
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));         // Check path
    }

    SUBCASE("Correct MSA with an deletion at the beginning") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "CATACCCGTAAAAGGG";
        std::vector<char> expected_cigar = {'D','D','D','M','M','M','M','M','M','M','M','M','M','M','M','M'};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1,false); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, 1, false);

        // Check if alignment was successful
        CHECK(alignment.compute_affine_gap_score(penalties) == 6); // Check score
        CHECK(alignment.edit_op == expected_cigar);                  // Check CIGAR
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));        // Check path
    }

    SUBCASE("Correct MSA with an insertion in the middle") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "ACCCGAAGGG";
        std::vector<char> expected_cigar = {'M','M','M','M','M','I','I','I','M','M','M','M','M'};

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1,false); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, false);

        // Check if alignment was successful
        CHECK(alignment.compute_affine_gap_score(penalties) == 6); // Check score
        CHECK(alignment.edit_op == expected_cigar);                  // Check CIGAR
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));        // Check path
    }

    SUBCASE("Correct MSA with diverging sequence") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string new_seq = "ACCCCCATAAGAGGG";

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1,false); // Create aligner

        // Align the new sequence
        theseus::Alignment alignment = aligner.align(new_seq, 1, false);

        // Check if alignment was successful
        CHECK(alignment.compute_affine_gap_score(penalties) == 9); // Check score
        CHECK(alignment.path == std::vector<NodeId>({0, 1, 2}));        // Check path
    }

    SUBCASE("Correct MSA with several sequences") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string seq_1       = "ACCCGTCAAAGGG";
        std::string seq_2       = "ACCCGAAGGG";
        std::string seq_3       = "ACCCGTCAAAGGG";
        std::string seq_4       = "ACCCCCATAAGAGGG";

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1,false); // Create aligner

        // Align and check sequence 1
        theseus::Alignment alignment = aligner.align(seq_1, false);
        CHECK(alignment.compute_affine_gap_score(penalties) == 2); // Check score

        // Align and check sequence 2
        alignment = aligner.align(seq_2, false);
        CHECK(alignment.compute_affine_gap_score(penalties) == 6); // Check score

        // Align and check sequence 3
        alignment = aligner.align(seq_3, false);
        CHECK(alignment.compute_affine_gap_score(penalties) == 0); // Check score

        // Align and check sequence 4
        alignment = aligner.align(seq_4, false);
        CHECK(alignment.compute_affine_gap_score(penalties) == 9); // Check score
    }

    SUBCASE("Multi-segment constructor and align_from") {
        // Build a 3-segment backbone graph: seg0="AAAA", seg1="CCCC", seg2="GGGG"
        std::vector<std::string_view> segments = {"AAAA", "CCCC", "GGGG"};
        std::vector<NodeId> node_ids;

        theseus::Penalties penalties(0, 2, 3, 1);
        theseus::Heuristics heuristics;
        theseus::TheseusMSA aligner(penalties, heuristics, segments, node_ids, 1);

        // node_ids should have 3 entries (one per segment)
        CHECK(node_ids.size() == 3);

        // Full alignment from source: sequence covers all 3 segments exactly
        theseus::Alignment alignment = aligner.align("AAAACCCCGGGG", 1, false, false);
        CHECK(alignment.compute_affine_gap_score(penalties) == 0);

        // align_from at seg1: sequence covers seg1+seg2 exactly
        alignment = aligner.align_from("CCCCGGGG", node_ids[1], 1, true);
        CHECK(alignment.compute_affine_gap_score(penalties) == 0);

        // align_from at seg0 with a mismatch in seg1
        alignment = aligner.align_from("AAAATTTTGGGG", node_ids[0], 1, false);
        CHECK(alignment.compute_affine_gap_score(penalties) == 8); // 4 mismatches * 2
    }

    SUBCASE("align_from partial coverage (sequence shorter than remaining graph)") {
        // Build a 4-segment backbone: seg0="AAAA", seg1="CCCC", seg2="GGGG", seg3="TTTT"
        std::vector<std::string_view> segments = {"AAAA", "CCCC", "GGGG", "TTTT"};
        std::vector<NodeId> node_ids;

        theseus::Penalties penalties(0, 2, 3, 1);
        theseus::Heuristics heuristics;
        theseus::TheseusMSA aligner(penalties, heuristics, segments, node_ids, 1);

        CHECK(node_ids.size() == 4);

        // Align from seg1 with a sequence that covers only seg1+seg2 (not seg3).
        // The aligner should handle partial coverage via init_partial_backtrace.
        theseus::Alignment alignment = aligner.align_from("CCCCGGGG", node_ids[1], 1, true);

        // Count matches in the CIGAR.
        int match_count = 0;
        for (char op : alignment.edit_op) {
            if (op == 'M') match_count++;
        }
        CHECK(match_count == 8);
    }

    SUBCASE("Correct reversed MSA") {
        // Equal sequences
        std::string initial_seq = "ACCCGTAAAAGGG";
        std::string seq_1       = "ACCCGTCAAAGGG";
        std::string seq_2       = "ACCCGAAGGG";
        std::string seq_3       = "GTAAAAGGG";
        // std::string seq_3       = "ACCCGTCAAAGGG";
        // std::string seq_4       = "ACCCCCATAAGAGGG";

        // Set aligner's parameters
        theseus::Penalties penalties(0, 2, 3, 1);                        // Create penalties object
        theseus::Heuristics heuristics;                  // Create heuristics object
        theseus::TheseusMSA aligner(penalties, heuristics, initial_seq,1,false); // Create aligner

        // Align and check sequence 1
        theseus::Alignment alignment = aligner.align(seq_1, 1, true);
        CHECK(alignment.compute_affine_gap_score(penalties) == 2); // Check score

        // Align and check sequence 2
        alignment = aligner.align(seq_2, 1, false);
        CHECK(alignment.compute_affine_gap_score(penalties) == 6); // Check score

        // // Align and check sequence 3
        alignment = aligner.align(seq_3, 1, true, true);
        CHECK(alignment.compute_affine_gap_score(penalties) == 0); // Check score

        // // Align and check sequence 4
        // alignment = aligner.align(seq_4, false);
        // CHECK(alignment.compute_affine_gap_score(penalties) == 9); // Check score
    }

    SUBCASE("End-to-end terminates (max-score safety cap)") {
        // Regression: an end-to-end align_from whose query cannot reach the
        // exact end node previously looped forever (no max-score cap with
        // heuristics off). It must now terminate (completed or unreachable),
        // never hang.
        std::vector<std::string_view> segments = {"AAAA", "CCCC", "GGGG", "TTTT"};
        std::vector<NodeId> node_ids;
        theseus::Penalties penalties(0, 2, 3, 1);
        theseus::Heuristics heuristics;
        theseus::TheseusMSA aligner(penalties, heuristics, segments, node_ids, 1);

        // End-to-end (is_ends_free=false), but the query is unrelated/short so
        // the exact end-node terminal is hard to hit. Must return, not hang.
        theseus::Alignment alignment =
            aligner.align_from("ACGTACGTACGT", node_ids[0], 1, /*is_ends_free=*/false,
                               /*start_offset=*/0, /*end_node=*/static_cast<int>(node_ids[2]));
        // No assertion on score; the point is that the call returns at all.
        CHECK(alignment.edit_op.size() >= 0);
    }

    SUBCASE("Deletion-heavy alignment does not read past query") {
        // Regression: the POA add path read new_seq[i] for X/D ops without
        // bounds-checking i against new_seq.size(), injecting garbage bytes.
        // A short query against a long backbone forces many deletions; under
        // ASan this would trip on the out-of-bounds read before the fix.
        std::vector<std::string_view> segments = {"AAAA", "CCCC", "GGGG", "TTTT"};
        std::vector<NodeId> node_ids;
        theseus::Penalties penalties(0, 2, 3, 1);
        theseus::Heuristics heuristics;
        theseus::TheseusMSA aligner(penalties, heuristics, segments, node_ids, 1);

        // Very short query spanning from seg0 toward the end: lots of deletions.
        theseus::Alignment alignment =
            aligner.align_from("AC", node_ids[0], 1, /*is_ends_free=*/true);
        // Every emitted base must be a valid nucleotide or gap (no garbage).
        for (char op : alignment.edit_op) {
            CHECK((op == 'M' || op == 'X' || op == 'I' || op == 'D'));
        }
    }
}