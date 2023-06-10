/**
  * FileName: node.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date:
  * Description: octree node
***/

#ifndef POISSON_RECONSTRUCTION_NODE_H
#define POISSON_RECONSTRUCTION_NODE_H

#include <iostream>
#include "mesh.h"

namespace poisson_reconstruction {

	// Octree node structure
	// With memory in mind, width and center are not stored in node
	// Dynamically allocates memory for the eight child nodes
	class Node {
	public:
		static const double epsilon;
		// Specifies the depth at which the cube is located
		uint8_t depth{ 0 };
		// Specifies the position of this cube node relative to the octree cube
		struct Offset {
			uint16_t x, y, z;
		}offset{ 0, 0, 0 };
		// Eight children Nodes
		Node* children{ nullptr };
		Node* parent{ nullptr };

		// when estimating local sampling density needs
		double weight_alpha{ 0.0 };

		// when computing the vector field needs
		bool set_vector_field{ false };
		double vector_field_factor[3]{ 0.0, 0.0, 0.0 };

		// the projection of the divergence of the vector field
		double proj_div_vec_field{ 0.0 };

		// used when compute the projection of indicator function
		uint32_t index{ 0 };

		// The solution component to a system of linear equations
		double solution_component{ 0.0 };

		// The type of intersection between the cube and the contour surface
		uint8_t cube_index;
		double corner_values[8];

		Node() = default;
		~Node() { if (children) delete[] children; }

		void init_eight_children();

		Position<double> get_center() const;
		double get_width() const;
		uint8_t get_depth() const { return depth; }

		// Determine whether two nodes intersect within a specified range
		static bool intersect(const Node* node1, double radius_coefficient1,
			const Node* node2, double radius_coefficient2);
		// Determine which of a parent_node's child nodes intersect with another_node
		// return a bitmask
		static uint8_t children_intersect(const Node* parent_node, double radius_coefficient1,
			const Node* another_node, double radius_coefficient2);
		// For each node that meets the condition in the octree, 
		// the operation F is performed between node and node2  
		template <typename F>
		static void top_down_process_node_to_node(Node* root, double radius_coefficient1,
			Node* node2, double radius_coefficient2,
			F& f);

		// Determine whether the node covers point within a specified range 
		static bool cover(const Node* node, double radius_coefficient, const Position<double>& point);
		// For each node that meets the condition in the octree, 
		// the operation F is performed between node and point  
		template <typename F>
		static void top_down_process_node_to_point(Node* root, double radius_coefficient,
			const Position<double>& point, F& f);

		// Gets the node adjacent to the parameter node's face
		static Node* node_adjacent_to_face(const Node* node, uint8_t face);
		// Gets the node adjacent to the parameter node's face
		static Node* node_adjacent_to_edge(const Node* node, uint8_t edge);

		static void faces_adjacent_to_edge(uint8_t edge, uint8_t& face1, uint8_t& face2);
		static void corner_adjacent_to_edge(uint8_t edge, uint8_t& corner1, uint8_t& corner2);
		static void corner_adjacent_to_face(uint8_t face, uint8_t& c1, uint8_t& c2, uint8_t& c3, uint8_t& c4);

		// clip empty leaf nodes
		static bool clip_tree(Node* root);

		Node* get_next_node(Node* node = nullptr);
		Node* get_next_branch_node(Node* node);
		Node* get_next_leaf_node(Node* node = nullptr);

		static long long get_corner_index(const Node* node, uint8_t corner, uint16_t index[3], uint8_t max_depth);
		static void get_corner_coordinate(const Node* node, uint8_t corner, Position<double>& point);
		static long long get_center_index(const Node* node, uint16_t index[3], uint8_t max_depth);

		static int total_nodes_count(const Node* root);
		static int total_leaf_node_count(const Node* root);
	};

	// Gets the center coordinates of the cube
	inline Position<double> Node::get_center() const {
		Position<double> ans;
		double w = get_width();
		ans.x = (offset.x + 0.5) * w;
		ans.y = (offset.y + 0.5) * w;
		ans.z = (offset.z + 0.5) * w;
		return ans;
	}

	inline double Node::get_width() const {
		return 1.0 / (1 << depth);
	}

	// assume that node1 and node2 are not nullptr
	// The radiation range of node1 is node1.width*radius_coefficient1
	// The radiation range of node2 is node2.width*radius_coefficient2
	inline bool Node::intersect(const Node* node1, double radius_coefficient1,
		const Node* node2, double radius_coefficient2)
	{
		double dist = node1->get_width() * radius_coefficient1 + node2->get_width() * radius_coefficient2;
		Position<double> center1 = node1->get_center();
		Position<double> center2 = node2->get_center();
		return fabs(center1.x - center2.x) < dist && fabs(center1.y - center2.y) < dist && fabs(center1.z - center2.z) < dist;
	}

	// assume that parent_node and another_node are not nullptr
	// and parent_node has eight children
	inline uint8_t Node::children_intersect(const Node* parent_node, double radius_coefficient1,
		const Node* another_node, double radius_coefficient2)
	{
		uint8_t bitmask = 0;
		for (uint8_t i = 0; i < 8; ++i) {
			if (intersect(&parent_node->children[i], radius_coefficient1, another_node, radius_coefficient2))
				bitmask |= (1 << i);
		}
		return bitmask;
	}

	inline bool Node::cover(const Node* node, double radius_coefficient, const Position<double>& point)
	{
		double dist = node->get_width() * radius_coefficient;
		Position<double> center = node->get_center();
		return fabs(center.x - point.x) < dist && fabs(center.y - point.y) < dist && fabs(center.z - point.z) < dist;
	}

	inline void Node::faces_adjacent_to_edge(uint8_t edge, uint8_t& face1, uint8_t& face2) {
		static const uint8_t adjacent_faces[12][2] = {
			{2, 4},  // edge 0
			{3, 4},  // edge 1
			{2, 5},  // edge 2
			{3, 5},  // edge 3
			{0, 4},  // edge 4
			{1, 4},  // edge 5
			{0, 5},  // edge 6
			{1, 5},  // edge 7
			{0, 2},  // edge 8
			{1, 2},  // edge 9
			{0, 3},  // edge 10
			{1, 3},  // edge 11
		};
		face1 = adjacent_faces[edge][0];
		face2 = adjacent_faces[edge][1];
	}

	inline void Node::corner_adjacent_to_edge(uint8_t edge, uint8_t& corner1, uint8_t& corner2) {
		static const uint8_t adjacent_corners[12][2] = {
			{0, 1},  // edge 0
			{2, 3},  // edge 1
			{4, 5},  // edge 2
			{6, 7},  // edge 3
			{0, 2},  // edge 4
			{1, 3},  // edge 5
			{4, 6},  // edge 6
			{5, 7},  // edge 7
			{0, 4},  // edge 8
			{1, 5},  // edge 9
			{2, 6},  // edge 10
			{3, 7},  // edge 11
		};
		corner1 = adjacent_corners[edge][0];
		corner2 = adjacent_corners[edge][1];
	}

	inline void Node::corner_adjacent_to_face(uint8_t face, uint8_t& c1, uint8_t& c2, uint8_t& c3, uint8_t& c4) {
		static const uint8_t corners_adjacent_face[6][4] = {
			{4, 6, 0, 2},  // face 0
			{5, 7, 1, 3},  // face 1
			{4, 5, 0, 1},  // face 2
			{6, 7, 2, 3},  // face 3
			{2, 3, 0, 1},  // face 4
			{4, 5, 6, 7}   // face 5
		};
		c1 = corners_adjacent_face[face][0];
		c2 = corners_adjacent_face[face][1];
		c3 = corners_adjacent_face[face][2];
		c4 = corners_adjacent_face[face][3];
	}

	// Starting from the root node, for every node that satisfies the given condition,
	// perform the function operation F between node and node2
	// assume that root and node2 are not nullptr
	template<typename F>
	void Node::top_down_process_node_to_node(Node* root, double radius_coefficient1,
		Node* node2, double radius_coefficient2,
		F& f)
	{
		if (intersect(root, radius_coefficient1, node2, radius_coefficient2)) f(root, node2);
		else return;
		if (root->children == nullptr) return;
		else {
			for (uint8_t i = 0; i < 8; ++i) {
				top_down_process_node_to_node(&root->children[i], radius_coefficient1, node2, radius_coefficient2, f);
			}
		}
	}

	// Starting from the root node, for every node that satisfies
	// the given condition: the node's scope covers the point,
	// perform the function operation F between node and point
	// assume that root is not nullptr
	template<typename F>
	void Node::top_down_process_node_to_point(Node* root, double radius_coefficient,
		const Position<double>& point, F& f)
	{
		if (cover(root, radius_coefficient, point)) f(root, point);
		else return;
		if (root->children == nullptr) return;
		else {
			for (uint8_t i = 0; i < 8; ++i) {
				top_down_process_node_to_point(&root->children[i], radius_coefficient, point, f);
			}
		}
	}

	inline long long Node::get_corner_index(const Node* node, uint8_t corner, uint16_t index[3], uint8_t max_depth) {
		uint16_t offset[3] = { node->offset.x, node->offset.y, node->offset.z };
		uint16_t pos[3] = { corner & 1, (corner & 2) >> 1, (corner & 4) >> 2 };
		for (uint8_t dir = 0; dir < 3; ++dir) {
			index[dir] = (offset[dir] + pos[dir]) << (max_depth + 1 - node->depth);
		}
		return (long long)index[0] | ((long long)index[1] << 15) | ((long long)index[2] << 30);
	}

	inline void Node::get_corner_coordinate(const Node* node, uint8_t corner, Position<double>& point) {
		uint16_t pos[3] = { corner & 1, (corner & 2) >> 1, (corner & 4) >> 2 };
		auto center = node->get_center();
		double temp = node->get_width() / 2.0;
		if (pos[0]) point.x = center.x + temp;
		else point.x = center.x - temp;
		if (pos[1]) point.y = center.y + temp;
		else point.y = center.y - temp;
		if (pos[2]) point.z = center.z + temp;
		else point.z = center.z - temp;
	}

	inline long long Node::get_center_index(const Node* node, uint16_t index[3], uint8_t max_depth) {
		uint16_t offset[3] = { node->offset.x, node->offset.y, node->offset.z };
		for (uint8_t dir = 0; dir < 3; ++dir) {
			index[dir] = ((offset[dir] << 1) + 1) << (max_depth + 1 - (node->depth + 1));
		}
		return (long long)index[0] | ((long long)index[1] << 15) | ((long long)index[2] << 30);
	}

}

#endif // !POISSON_RECONSTRUCTION_NODE_H
