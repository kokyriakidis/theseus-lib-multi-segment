#include "../include/theseus/graph.h"

#include <algorithm>
#include <string>
#include <vector>
#include <ranges>

namespace theseus {

class GraphImpl {
public:
    using NodeId = Graph::NodeId;
    using NodeView = Graph::NodeView;
    using NodeIdRange = Graph::NodeIdRange;

    /*
     * Forward iterator for alive nodes in the graph. This iterator skips over
     * nodes that have been removed (i.e., marked as not alive) from the graph.
     */
    class AliveNodeIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = NodeId;
        using difference_type = std::ptrdiff_t;
        using pointer = const NodeId *;
        using reference = const NodeId &;

        AliveNodeIterator(const GraphImpl *graph, NodeId current)
            : graph_(graph), current_(current) {
            advance_to_alive();
        }

        reference operator*() const { return current_; }
        pointer operator->() const { return &current_; }

        AliveNodeIterator &operator++() {
            ++current_;
            advance_to_alive();
            return *this;
        }

        AliveNodeIterator operator++(int) {
            AliveNodeIterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const AliveNodeIterator &other) const {
            return graph_ == other.graph_ && current_ == other.current_;
        }

    private:
        void advance_to_alive() {
            while (current_ < graph_->nodes_.size() && !graph_->nodes_[current_].alive) {
                ++current_;
            }
        }

        const GraphImpl *graph_;
        NodeId current_;
    };

    /**
     * Copy the contents of another GraphImpl into this one. Explicit function
     * instead of copy constructor/assignment operator used by Graph.
     *
     * @param other The GraphImpl to copy from.
     */
    void copy_from(const GraphImpl& other) {
        nodes_ = other.nodes_;
        free_node_ids_ = other.free_node_ids_;
        source_nodes_ = other.source_nodes_;
        sink_nodes_ = other.sink_nodes_;
    }

    NodeId add_node(std::string_view sequence) {
        return add_node(std::string(sequence));
    }

    NodeId add_node(std::string&& sequence) {
        NodeId id = 0;

        if (!free_node_ids_.empty()) {
            id = free_node_ids_.back();
            free_node_ids_.pop_back();

            nodes_[id].sequence = std::move(sequence);
            nodes_[id].alive = true;
        }
        else {
            id = nodes_.size();

            nodes_.emplace_back(Node{std::move(sequence), {}, {}, true});
        }

        source_nodes_.push_back(id);
        sink_nodes_.push_back(id);

        return id;
    }

    void expand_sequence(NodeId id, std::string_view suffix) {
        if (!is_valid_node(id)) {
            throw Graph::InvalidNodeException(id);
        }

        nodes_[id].sequence.append(suffix);
    }

    std::string split_sequence(NodeId id, size_t idx) {
        if (!is_valid_node(id)) {
            throw Graph::InvalidNodeException(id);
        }

        std::string new_seq = nodes_[id].sequence.substr(idx);
        nodes_[id].sequence.resize(idx);

        return new_seq;
    }

    void remove_node(NodeId id) {
        if (!is_valid_node(id)) {
            throw Graph::InvalidNodeException(id);
        }

        if (is_source(id)) {
            remove_id_from_vec(source_nodes_, id);
        }

        if (is_sink(id)) {
            remove_id_from_vec(sink_nodes_, id);
        }

        for (NodeId out_id : nodes_[id].out_nodes) {
            remove_id_from_vec(nodes_[out_id].in_nodes, id);
        }

        for (NodeId in_id : nodes_[id].in_nodes) {
            remove_id_from_vec(nodes_[in_id].out_nodes, id);
        }

        nodes_[id].sequence.clear();
        nodes_[id].in_nodes.clear();
        nodes_[id].out_nodes.clear();
        nodes_[id].alive = false;

        free_node_ids_.push_back(id);
    }

    void add_edge(NodeId from, NodeId to) {
        if (!is_valid_node(from)) {
            throw Graph::InvalidNodeException(from);
        }
        if (!is_valid_node(to)) {
            throw Graph::InvalidNodeException(to);
        }

        nodes_[from].out_nodes.push_back(to);
        nodes_[to].in_nodes.push_back(from);

        remove_id_from_vec(source_nodes_, to);
        remove_id_from_vec(sink_nodes_, from);
    }

    bool remove_edge(NodeId from, NodeId to) {
        if (!is_valid_node(from)) {
            throw Graph::InvalidNodeException(from);
        }
        if (!is_valid_node(to)) {
            throw Graph::InvalidNodeException(to);
        }

        bool removed = remove_id_from_vec(nodes_[from].out_nodes, to);
        if (removed) {
            remove_id_from_vec(nodes_[to].in_nodes, from);

            if (nodes_[from].out_nodes.empty()) {
                sink_nodes_.push_back(from);
            }
            if (nodes_[to].in_nodes.empty()) {
                source_nodes_.push_back(to);
            }
        }

        return removed;
    }

    void remove_out_edges(NodeId from) {
        if (!is_valid_node(from)) {
            throw Graph::InvalidNodeException(from);
        }

        for (NodeId to : nodes_[from].out_nodes) {
            remove_id_from_vec(nodes_[to].in_nodes, from);

            if (nodes_[to].in_nodes.empty()) {
                source_nodes_.push_back(to);
            }
        }

        nodes_[from].out_nodes.clear();

        if (!is_sink(from)) {
            sink_nodes_.push_back(from);
        }
    }

    void remove_in_edges(NodeId to) {
        if (!is_valid_node(to)) {
            throw Graph::InvalidNodeException(to);
        }

        for (NodeId from : nodes_[to].in_nodes) {
            remove_id_from_vec(nodes_[from].out_nodes, to);

            if (nodes_[from].out_nodes.empty()) {
                sink_nodes_.push_back(from);
            }
        }

        nodes_[to].in_nodes.clear();

        if (!is_source(to)) {
            source_nodes_.push_back(to);
        }
    }

    void substitute_out_edges(NodeId orig_from, NodeId new_from) {
        if (!is_valid_node(orig_from)) {
            throw Graph::InvalidNodeException(orig_from);
        }
        if (!is_valid_node(new_from)) {
            throw Graph::InvalidNodeException(new_from);
        }

        for (NodeId to : nodes_[orig_from].out_nodes) {
            remove_id_from_vec(nodes_[to].in_nodes, orig_from);
            nodes_[to].in_nodes.push_back(new_from);
        }

        nodes_[new_from].out_nodes = std::move(nodes_[orig_from].out_nodes);
        nodes_[orig_from].out_nodes.clear();

        if (!is_sink(orig_from)) {
            sink_nodes_.push_back(orig_from);
        }
        if (!nodes_[new_from].out_nodes.empty()) {
            remove_id_from_vec(sink_nodes_, new_from);
        }
    }

    void substitute_in_edges(NodeId orig_to, NodeId new_to) {
        if (!is_valid_node(orig_to)) {
            throw Graph::InvalidNodeException(orig_to);
        }
        if (!is_valid_node(new_to)) {
            throw Graph::InvalidNodeException(new_to);
        }

        for (NodeId from : nodes_[orig_to].in_nodes) {
            remove_id_from_vec(nodes_[from].out_nodes, orig_to);
            nodes_[from].out_nodes.push_back(new_to);
        }

        nodes_[new_to].in_nodes = std::move(nodes_[orig_to].in_nodes);
        nodes_[orig_to].in_nodes.clear();

        if (!is_source(orig_to)) {
            source_nodes_.push_back(orig_to);
        }
        if (!nodes_[new_to].in_nodes.empty()) {
            remove_id_from_vec(source_nodes_, new_to);
        }
    }

    NodeView node(NodeId id) const {
        if (!is_valid_node(id)) {
            throw Graph::InvalidNodeException(id);
        }

        const auto &node = nodes_[id];
        return {
            Graph::SequenceView(node.sequence, false),
            Graph::NodeIdRange(node.in_nodes.cbegin(), node.in_nodes.cend()),
            Graph::NodeIdRange(node.out_nodes.cbegin(), node.out_nodes.cend()),
        };
    }

    NodeView node_rev(NodeId id) const {
        if (!is_valid_node(id)) {
            throw Graph::InvalidNodeException(id);
        }

        const auto &node = nodes_[id];
        return {
            Graph::SequenceView(node.sequence, true),
            Graph::NodeIdRange(node.out_nodes.cbegin(), node.out_nodes.cend()),
            Graph::NodeIdRange(node.in_nodes.cbegin(), node.in_nodes.cend()),
        };
    }

    size_t nnodes() const {
        return nodes_.size() - free_node_ids_.size();
    }

    size_t node_id_bound() const {
        return nodes_.size();
    }

    NodeIdRange nodes() const {
        return NodeIdRange(AliveNodeIterator(this, 0),
                           AliveNodeIterator(this, nodes_.size()));
    }

    bool is_source(NodeId id) const {
        if (!is_valid_node(id)) {
            throw Graph::InvalidNodeException(id);
        }

        return nodes_[id].in_nodes.empty();
    }

    bool is_sink(NodeId id) const {
        if (!is_valid_node(id)) {
            throw Graph::InvalidNodeException(id);
        }
        return nodes_[id].out_nodes.empty();
    }

    int node_size(NodeId id) const {
        if (!is_valid_node(id)) {
            throw Graph::InvalidNodeException(id);
        }
        return nodes_[id].sequence.size();
    }

    NodeIdRange source_nodes() const {
        return NodeIdRange(source_nodes_.cbegin(), source_nodes_.cend());
    }

    NodeIdRange sink_nodes() const {
        return NodeIdRange(sink_nodes_.cbegin(), sink_nodes_.cend());
    }

    void print_graph() {
        std::cout << "\\begin{tikzpicture}[\n";
        std::cout << "  >=stealth,\n";
        std::cout << "  every node/.style={draw, rounded corners, minimum height=8mm, font=\\ttfamily\\small}\n";
        std::cout << "]\n";

        for (long unsigned int i = 0; i < nodes_.size(); ++i) {
            auto node = nodes_[i];
            std::cout << "  \\node (n" << i << ") at (" << (i * 3) << ",0) {" << i << ": " << node.sequence << "};\n";
        }

        for (long unsigned int i = 0; i < nodes_.size(); ++i) {
            auto node = nodes_[i];
            for (auto out_id : node.out_nodes) {
                std::cout << "  \\draw[->] (n" << i << ") -- (n" << out_id << ");\n";
            }
        }

        std::cout << "\\end{tikzpicture}\n";
    }

private:
    struct Node {
        std::string sequence;
        std::vector<NodeId> in_nodes;
        std::vector<NodeId> out_nodes;
        bool alive;
    };

    std::vector<Node> nodes_;

    std::vector<NodeId> free_node_ids_;

    std::vector<NodeId> source_nodes_;
    std::vector<NodeId> sink_nodes_;

    /**
     * Remove the given id from the vector if it exists. Avoid dynamic allocations
     * by swapping the id to remove with the last element and popping the back.
     *
     * @param vec The vector from which to remove the id.
     * @param id The id to remove.
     * @return true if the id was found and removed, false otherwise.
     */
    static bool remove_id_from_vec(std::vector<NodeId> &vec, NodeId id) {
        auto it = std::find(vec.begin(), vec.end(), id);
        if (it != vec.end()) {
            *it = vec.back();
            vec.pop_back();
            return true;
        }
        return false;
    }

    /**
     * Check if the given NodeId corresponds to a valid node in the graph.
     *
     * @param id The NodeId to check.
     * @return true if the NodeId is valid, false otherwise.
     */
    bool is_valid_node(NodeId id) const {
        return id < nodes_.size() && nodes_[id].alive;
    }
};

Graph::Graph() : impl_(std::make_unique<GraphImpl>()) {}

Graph::~Graph() = default;

Graph::Graph(const Graph& other) : impl_(std::make_unique<GraphImpl>()) {
    impl_->copy_from(*other.impl_);
}

Graph& Graph::operator=(const Graph& other) {
    if (this != &other) {
        impl_->copy_from(*other.impl_);
    }
    return *this;
}

Graph::Graph(Graph&&) noexcept = default;

Graph& Graph::operator=(Graph&&) noexcept = default;

void Graph::swap(Graph& other) noexcept {
    impl_.swap(other.impl_);
}

Graph::NodeId Graph::add_node(std::string_view sequence) {
    return impl_->add_node(sequence);
}

Graph::NodeId Graph::add_node(const char *sequence) {
    return impl_->add_node(std::string_view(sequence));
}

Graph::NodeId Graph::add_node(std::string &&sequence) {
    return impl_->add_node(std::move(sequence));
}

void Graph::expand_sequence(Graph::NodeId id, std::string_view suffix) {
    impl_->expand_sequence(id, suffix);
}

std::string Graph::split_sequence(NodeId id, size_t idx) {
    return impl_->split_sequence(id, idx);
}

void Graph::remove_node(NodeId id) { impl_->remove_node(id); }

void Graph::add_edge(NodeId from, NodeId to) { impl_->add_edge(from, to); }

bool Graph::remove_edge(NodeId from, NodeId to) {
    return impl_->remove_edge(from, to);
}

void Graph::remove_out_edges(NodeId from) { impl_->remove_out_edges(from); }

void Graph::remove_in_edges(NodeId to) { impl_->remove_in_edges(to); }

void Graph::substitute_out_edges(NodeId orig_from, NodeId new_from) {
    impl_->substitute_out_edges(orig_from, new_from);
}

void Graph::substitute_in_edges(NodeId orig_to, NodeId new_to) {
    impl_->substitute_in_edges(orig_to, new_to);
}

Graph::NodeView Graph::node(NodeId id) const { return impl_->node(id); }

Graph::NodeView Graph::node_rev(NodeId id) const { return impl_->node_rev(id); }

size_t Graph::nnodes() const { return impl_->nnodes(); }
size_t Graph::node_id_bound() const { return impl_->node_id_bound(); }

bool Graph::is_source(NodeId id) const { return impl_->is_source(id); }

bool Graph::is_sink(NodeId id) const { return impl_->is_sink(id); }

Graph::NodeIdRange Graph::nodes() const { return impl_->nodes(); }

int Graph::node_size(NodeId id) const { return impl_->node_size(id); }

Graph::NodeIdRange Graph::source_nodes() const { return impl_->source_nodes(); }

Graph::NodeIdRange Graph::sink_nodes() const { return impl_->sink_nodes(); }

void Graph::print_graph() { return impl_->print_graph(); }

}  // namespace theseus