/**
  * FileName: node.cpp
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date:
  * Description: octree node
***/

#include "node.h"

namespace poisson_reconstruction {

	//---------------------------------------------------------------------------------------
	// Node member functions
	//---------------------------------------------------------------------------------------

	const double Node::epsilon = 1e-6;

	void Node::init_eight_children() {
		if (children) {
			delete[] children;
			children = nullptr;
		}

		Offset child_offset;
		uint8_t child_number;
		children = new Node[8];
		for (uint8_t x = 0; x < 2; ++x)
			for (uint8_t y = 0; y < 2; ++y)
				for (uint8_t z = 0; z < 2; ++z) {
					// handle the child with number zyx, i.e at (x, y, z)
					child_number = (z << 2) | (y << 1) | x;
					// set children's Offset
					child_offset.x = (offset.x << 1) + x;
					child_offset.y = (offset.y << 1) + y;
					child_offset.z = (offset.z << 1) + z;
					children[child_number].offset = child_offset;
					children[child_number].depth = depth + 1;
					children[child_number].parent = this;
				}
	}

	// assume that node is not nullptr
	Node* Node::node_adjacent_to_face(const Node* node, uint8_t face)
	{	
		if (node->parent == nullptr) return nullptr;

		// Check whether the neighboring(answer) node is the brother of this node
		// 3 internal faces of this node relative to the parent node:
		// this node's number = node - parent->children (parent = node->parent)
		// number <-> (x, y, z)
		// x-direction: 00(x^1)
		// y-direction: 01(y^1)
		// z-direction: 10(z^1)
		uint8_t dir = (face & 6) >> 1; // 6 <-> 110 
		const Node* parent = node->parent;
		uint8_t number = node - parent->children;  // zyx

		// compute the number of the neighboring(answer) node 
		Position<uint8_t> coor(number & 1, (number & 2) >> 1, (number & 4) >> 2);
		uint8_t ans_number;
		switch (dir) {
		case 0:  // x-direction
			ans_number = coor.z << 2 | coor.y << 1 | (coor.x ^ 1);
			break;
		case 1:  // y-direction 
			ans_number = coor.z << 2 | (coor.y ^ 1) << 1 | coor.x;
			break;
		case 2:  // z-direction
			ans_number = (coor.z ^ 1) << 2 | coor.y << 1 | coor.x;
			break;
		}

		if (face == (0 | (coor.x ^ 1)) || face == (2 | (coor.y ^ 1)) || face == (4 | (coor.z ^ 1))) {
			return &parent->children[ans_number];
		}

		// Search recursively to find the node adjacent to the parent node on this parameter face
		Node* parent_adjacent_node = node_adjacent_to_face(parent, face);
		if (parent_adjacent_node == nullptr || parent_adjacent_node->children == nullptr) return nullptr;

		return &parent_adjacent_node->children[ans_number];
	}

	Node* Node::node_adjacent_to_edge(const Node* node, uint8_t edge)
	{
		if (node->parent == nullptr) return nullptr;

		Node* parent = node->parent;
		uint8_t dir = (edge & 12) >> 2;  // 12 <-> 1100
		uint8_t number = node - parent->children;

		// compute the number of the neighboring(answer) node 
		Position<uint8_t> coor(number & 1, (number & 2) >> 1, (number & 4) >> 2);
		uint8_t ans_number;
		switch (dir) {
		case 0:  // x-direction
			ans_number = (coor.z ^ 1) << 2 | (coor.y ^ 1) << 1 | coor.x;
			break;
		case 1:  // y-direction
			ans_number = (coor.z ^ 1) << 2 | coor.y << 1 | (coor.x ^ 1);
			break;
		case 2:  // z-direction
			ans_number = coor.z << 2 | (coor.y ^ 1) << 1 | (coor.x ^ 1);
			break;
		}

		// these edges are inside the parent node cube
		// It's different for each of these 8 nodes
		static const uint8_t internal_edges[8][3] = {
			{3, 7, 11},  // child 0
			{3, 6, 10},  // child 1
			{2, 7,  9},  // child 2
			{2, 6,  8},  // child 3
			{1, 5, 11},  // child 4
			{1, 4, 10},  // child 5
			{0, 5,  9},  // child 6
			{0, 4,  8},  // child 7
		};
		// 1-th case: Check whether the neighboring(answer) node is the brother of this node
		if (edge == internal_edges[number][0] || edge == internal_edges[number][1] || edge == internal_edges[number][2]) {
			return &parent->children[ans_number];
		}

		// these edges are outside the parent node cube
		// It's different for each of these 8 nodes
		static const uint8_t external_edges[8][3] = {
			{0, 4,  8},  // child 0
			{0, 5,  9},  // child 1
			{1, 4, 10},  // child 2
			{1, 5, 11},  // child 3
			{2, 6,  8},  // child 4
			{2, 7,  9},  // child 5
			{3, 6, 10},  // child 6
			{3, 7, 11},  // child 7
		};
		// 2-th Check whether the neighboring(answer) node is the brother of this node
		if (edge == external_edges[number][0] || edge == external_edges[number][1] || edge == external_edges[number][2]) {
			// Search recursively to find the node adjacent to the parent node on this parameter edge
			Node* parent_adjacent_node = node_adjacent_to_edge(parent, edge);
			if (parent_adjacent_node == nullptr || parent_adjacent_node->children == nullptr) return nullptr;

			return &parent_adjacent_node->children[ans_number];
		}

		// these edges are in the surface of parent node cube
		// It's different for each of these 8 nodes
		static const uint8_t surface_edges[8][6] = {
			// x      y       z
			{6, 10,  2, 9,  1, 5},  // child 0
			{7, 11,  2, 8,  1, 4},  // child 1
			{6,  8,  3,11,  0, 5},  // child 2
			{7,  9,  3,10,  4, 0},  // child 3
			{4, 10,  0, 9,  3, 7},  // child 4
			{5, 11,  0, 8,  3, 6},  // child 5
			{4,  8,  1,11,  2, 7},  // child 6
			{5,  9,  1,10,  2, 6},  // child 7
		};
		for (uint8_t dir = 0; dir < 2; ++dir) {
			if (edge == surface_edges[number][2 * dir] || edge == surface_edges[number][2 * dir + 1]) {
				Node* parent_adjacent_node = nullptr;
				switch (dir) {
				case 0:  // face in x-direction
					// Search recursively to find the node adjacent to the parent node on specified face
					parent_adjacent_node = node_adjacent_to_face(parent, 0 | coor.x);
					break;
				case 1:  // face in y-direction
					parent_adjacent_node = node_adjacent_to_face(parent, 2 | coor.y);
					break;
				case 2:  // face in z-direction
					parent_adjacent_node = node_adjacent_to_face(parent, 4 | coor.z);
					break;
				}

				// 3-th case
				if (parent_adjacent_node == nullptr || parent_adjacent_node->children == nullptr) return nullptr;
				return &parent_adjacent_node->children[ans_number];
			}
		}

		// never reach here
		return nullptr;
	}

	// Determines whether subtree with parameter root as its root can be clipped
	bool Node::clip_tree(Node* root)
	{
		if (root == nullptr) return true;
		if (root->children == nullptr) return !root->set_vector_field;

		// clip from bottom to top
		bool clip = true;
		for (uint8_t i = 0; i < 8; ++i) {
			if (!clip_tree(&root->children[i])) {
				clip = false;
			}
		}
		// When all eight child nodes can be clipped
		if (clip && root->children != nullptr) {
			delete[] root->children;
			root->children = nullptr;
		}

		if (clip == false) return false;
		else return !root->set_vector_field;
	}

	Node* Node::get_next_node(Node* node) {
		// when call this func with root == nulllptr, return the caller node
		if (node == nullptr) return this;

		if (node->children != nullptr) {
			return &node->children[0];
		}

		return get_next_branch_node(node);
	}

	// get the next branch node
	// assume that node is not nullptr
	Node* Node::get_next_branch_node(Node* node) {
		if (node->parent == nullptr) return nullptr;

		Node* parent = node->parent;
		uint8_t node_number = node - parent->children;
		if (node_number < 7) return &parent->children[node_number + 1];

		// find the parent's next branch node
		return get_next_branch_node(parent);
	}

	Node* Node::get_next_leaf_node(Node* node) {
		Node* ans;
		// get the first leaf node
		if (node == nullptr) {
			ans = this;
			while (ans->children) ans = &ans->children[0];
			return ans;
		}
		ans = get_next_branch_node(node);
		if (ans == nullptr) return nullptr;
		while (ans->children) ans = &ans->children[0];
		return ans;
	}

	int Node::total_nodes_count(const Node* root) {
		if (root == nullptr) return 0;
		if (root->children == nullptr) return 1;

		int ans = 1;
		for (uint8_t i = 0; i < 8; ++i) {
			ans += total_nodes_count(&root->children[i]);
		}
		return ans;
	}

	int Node::total_leaf_node_count(const Node* root) {
		if (root == nullptr) return 0;
		if (root->children == nullptr) return 1;

		int ans = 0;
		for (uint8_t i = 0; i < 8; ++i) {
			ans += total_leaf_node_count(&root->children[i]);
		}
		return ans;
	}

}