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


#include <vector>
#include <string>

#include "graph.h"
#include "gfa_graph.h"

namespace theseus {
    // Constructor from GFA stream
    Graph::Graph(std::istream &gfa_stream) {
        GfaGraph gfa_graph(gfa_stream);

        // Add nodes
        _vertices.reserve(gfa_graph.gfa_nodes.size());
        for (int i = 0; i < gfa_graph.gfa_nodes.size(); ++i) {
            vertex v;
            v.name = gfa_graph.gfa_nodes[i].name;
            v.value = gfa_graph.gfa_nodes[i].seq;
            _vertices.push_back(v);
        }

        // Add edges
        for (int i = 0; i < gfa_graph.gfa_edges.size(); ++i) {
            edge e;
            e.from_vertex = gfa_graph.gfa_edges[i].from_node;
            e.to_vertex = gfa_graph.gfa_edges[i].to_node;
            e.overlap = gfa_graph.gfa_edges[i].overlap;
            _vertices[e.from_vertex].out_edges.push_back(e);
            _vertices[e.to_vertex].in_edges.push_back(e);
        }

        // Create name to id mapping
        for (int i = 0; i < _vertices.size(); ++i) {
            name_to_id_[_vertices[i].name] = i;
        }
    }

    // Constructor from handle graph
    Graph::Graph(const handlegraph::HandleGraph &handle_graph) {
        // Add nodes
        handle_graph.for_each_handle([&](const handlegraph::handle_t& h) {
            // Forward orientation
            vertex v;
            v.name = std::to_string(handle_graph.get_id(h)); // Use the handle id as the name of the vertex
            v.value = handle_graph.get_sequence(h);
            _vertices.push_back(v);

            // Reverse orientation
            const handlegraph::handle_t& h_rev = handle_graph.flip(h);
            v.name = std::to_string(handle_graph.get_id(h_rev));
            v.value = handle_graph.get_sequence(h_rev);
            _vertices.push_back(v);
        });

        // Create name to id mapping
        for (int i = 0; i < _vertices.size(); ++i) {
            name_to_id_[_vertices[i].name] = i;
        }

        // Add edges
        handle_graph.for_each_edge([&](const handlegraph::edge_t& e) {
            // Forward orientation
            edge edge_fwd;
            std::string from_handle_name = std::to_string(handle_graph.get_id(e.first));
            std::string to_handle_name   = std::to_string(handle_graph.get_id(e.second));
            int from_vertex_id = name_to_id_.at(from_handle_name);
            int to_vertex_id   = name_to_id_.at(to_handle_name);
            // Check that the vertices exist in the graph
            if (from_vertex_id == -1 || to_vertex_id == -1) {
                std::cerr << "Error: Vertex not found in the graph for edge (" << from_handle_name << " -> " << to_handle_name << ")" << std::endl;
                return;
            }
            edge_fwd.from_vertex = from_vertex_id;
            edge_fwd.to_vertex   = to_vertex_id;
            edge_fwd.overlap     = 0; // TODO: Now we suppose that overlap is 0
            _vertices[edge_fwd.from_vertex].out_edges.push_back(edge_fwd);
            _vertices[edge_fwd.to_vertex].in_edges.push_back(edge_fwd);

            // Reverse orientation
            edge edge_rev;
            std::string rev_from_handle_name = std::to_string(handle_graph.get_id(handle_graph.flip(e.second)));
            std::string rev_to_handle_name   = std::to_string(handle_graph.get_id(handle_graph.flip(e.first)));
            int rev_from_vertex_id = name_to_id_.at(rev_from_handle_name);
            int rev_to_vertex_id   = name_to_id_.at(rev_to_handle_name);
            // Check that the vertices exist in the graph
            if (rev_from_vertex_id == -1 || rev_to_vertex_id == -1) {
                std::cerr << "Error: Vertex not found in the graph for reverse edge (" << rev_from_handle_name << " -> " << rev_to_handle_name << ")" << std::endl;
                return;
            }
            edge_rev.from_vertex = rev_from_vertex_id;
            edge_rev.to_vertex   = rev_to_vertex_id;
            edge_rev.overlap     = 0; // TODO: Now we suppose that overlap is 0
            _vertices[edge_rev.from_vertex].out_edges.push_back(edge_rev);
            _vertices[edge_rev.to_vertex].in_edges.push_back(edge_rev);
        });
    }
}