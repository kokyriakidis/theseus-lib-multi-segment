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

namespace theseus {

TheseusMSA::TheseusMSA(const Penalties &penalties,
                       const Heuristics &heuristics,
                       std::string_view seq) {

    // Create the initial graph
    theseus::Graph G;
    theseus::Graph::vertex source_v, central_v, sink_v;
    Graph::edge source_edge, central_edge;

    // Source vertex
    source_edge.from_vertex = 0;
    source_edge.to_vertex = 1;
    source_v.out_edges.push_back(source_edge);
    source_v.first_poa_vtx = 0;
    G._vertices.push_back(source_v);

    // Central vertex (initial sequence)
    central_edge.from_vertex = 1;
    central_edge.to_vertex = 2;
    central_v.in_edges.push_back(source_edge);
    central_v.out_edges.push_back(central_edge);
    central_v.first_poa_vtx = 1;
    central_v.value = seq;
    G._vertices.push_back(central_v);

    // Sink vertex
    sink_v.in_edges.push_back(central_edge);
    sink_v.first_poa_vtx = seq.size() + 1;
    G._vertices.push_back(sink_v);

    msa_aligner_impl_ = std::make_unique<TheseusAlignerImpl>(penalties, heuristics, std::move(G), true);
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
 * @param start_node
 * @param start_offset
 * @return Alignment
 */
Alignment TheseusMSA::align(
    std::string_view seq) {

    std::string start_node;
    return msa_aligner_impl_->align(seq, start_node, 0);
}

/**
 * @brief Print the current POA graph in MSA format.
 *
 */
void TheseusMSA::print_as_gfa(std::ofstream &out_stream) {
    msa_aligner_impl_->print_as_gfa(out_stream);
}

/**
 * @brief Print the current POA graph in MSA format.
 *
 */
void TheseusMSA::print_as_msa(std::ofstream &out_stream) {
    msa_aligner_impl_->print_as_msa(out_stream);
}

/**
 * @brief Return consensus sequence.
 *
 */
std::string TheseusMSA::get_consensus_sequence() {
    return msa_aligner_impl_->get_consensus_sequence();
}

/**
 * @brief Print in graphviz format.
 *
 */
void TheseusMSA::print_as_dot(std::ofstream &out_stream) {
    msa_aligner_impl_->print_as_dot(out_stream);
}

} // namespace theseus