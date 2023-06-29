/**
  * FileName: magic_cube.cpp
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date:
  * Description:
***/

#include "magic_cube.h"
#include "mesh.h"

namespace poisson_reconstruction {

	//---------------------------------------------------------------------------------------
	// Octree::Magic_cube member functions
	//---------------------------------------------------------------------------------------

	Magic_cube::Magic_cube() {
		for (int x = 0; x < 2; ++x)
			for (int y = 0; y < 2; ++y)
				for (int z = 0; z < 2; ++z) {
					neighbor[x][y][z] = nullptr;
				}
	}

	void Magic_cube::construct(Node* cnode) {
		// The search direction of the parent node cube
		static const Position<int8_t> search_matrix[8][7] = {
			// for face   x            y                   z                for edge   x               y                  z                 for vertex
			{Position<int8_t>(-1, 0, 0), Position<int8_t>(0, -1, 0), Position<int8_t>(0, 0, -1),  Position<int8_t>(0, -1, -1), Position<int8_t>(-1, 0, -1), Position<int8_t>(-1, -1, 0),  Position<int8_t>(-1, -1, -1)}, // for child 0
			{Position<int8_t>(1, 0, 0), Position<int8_t>(0, -1, 0), Position<int8_t>(0, 0, -1),  Position<int8_t>(0, -1, -1), Position<int8_t>(1, 0, -1), Position<int8_t>(1,  1, 0),  Position<int8_t>(1, -1, -1)}, // for child 1
			{Position<int8_t>(-1, 0, 0), Position<int8_t>(0,  1, 0), Position<int8_t>(0, 0, -1),  Position<int8_t>(0,  1, -1), Position<int8_t>(-1, 0, -1), Position<int8_t>(-1,  1, 0),  Position<int8_t>(-1,  1, -1)}, // for child 2
			{Position<int8_t>(1, 0, 0), Position<int8_t>(0,  1, 0), Position<int8_t>(0, 0, -1),  Position<int8_t>(0,  1, -1), Position<int8_t>(1, 0, -1), Position<int8_t>(1,  1, 0),  Position<int8_t>(1,  1, -1)}, // for child 3
			{Position<int8_t>(-1, 0, 0), Position<int8_t>(0, -1, 0), Position<int8_t>(0, 0,  1),  Position<int8_t>(0, -1,  1), Position<int8_t>(-1, 0,  1), Position<int8_t>(-1, -1, 0),  Position<int8_t>(-1, -1,  1)}, // for child 4
			{Position<int8_t>(1, 0, 0), Position<int8_t>(0, -1, 0), Position<int8_t>(0, 0,  1),  Position<int8_t>(0, -1,  1), Position<int8_t>(1, 0,  1), Position<int8_t>(1, -1, 0),  Position<int8_t>(1, -1,  1)}, // for child 5
			{Position<int8_t>(-1, 0, 0), Position<int8_t>(0,  1, 0), Position<int8_t>(0, 0,  1),  Position<int8_t>(0,  1,  1), Position<int8_t>(-1, 0,  1), Position<int8_t>(-1,  1, 0),  Position<int8_t>(-1,  1,  1)}, // for child 6
			{Position<int8_t>(1, 0, 0), Position<int8_t>(0,  1, 0), Position<int8_t>(0, 0,  1),  Position<int8_t>(0,  1,  1), Position<int8_t>(1, 0,  1), Position<int8_t>(1,  1, 0),  Position<int8_t>(1,  1,  1)}, // for child 7
		};

		std::memset(neighbor, 0, sizeof(Node*) * 27);

		if (cnode->parent == nullptr) {
			neighbor[1][1][1] = cnode;
			return;
		}

		Node* parent = cnode->parent;
		// the number of cnode in its parent node
		uint8_t number = cnode - parent->children;  // zyx
		//                      x            y           z
		Position<int8_t> pos(number & 1, (number & 2) >> 1, (number & 4) >> 2);
		Position<int8_t> trans;

		// set 8 nodes of 27 in total
		// set cnode as the center at neighbor[1][1][1]
		trans = Position<int8_t>(pos.x ^ 1, pos.y ^ 1, pos.z ^ 1);
		for (uint8_t x = 0; x < 2; ++x)
			for (uint8_t y = 0; y < 2; ++y)
				for (uint8_t z = 0; z < 2; ++z) {
					neighbor[x + trans.x][y + trans.y][z + trans.z] = &parent->children[z << 2 | y << 1 | x];
				}

		Magic_cube parent_magic_cube;
		parent_magic_cube.construct(parent);

		Node* temp;
		Position<int8_t> face_adjacent_node_pos;
		Position<int8_t> center(1, 1, 1);
		// set 12 nodes of 27
		// Find 12 cube nodes adjacent to cnode's 3 faces
		// this 3 faces are: e.g  0 -> face 0, 2, 4;  ..., 5 -> face 1, 2, 5, ...  
		for (uint8_t dir = 0; dir < 3; ++dir) {
			temp = parent_magic_cube(center + search_matrix[number][dir]);
			if (temp != nullptr) {
				// when needed
				if (temp->children == nullptr) temp->init_eight_children();

				switch (dir) {
				case 0: face_adjacent_node_pos = Position<int8_t>(pos.x ^ 1, pos.y, pos.z); break;
				case 1: face_adjacent_node_pos = Position<int8_t>(pos.x, pos.y ^ 1, pos.z); break;
				case 2: face_adjacent_node_pos = Position<int8_t>(pos.x, pos.y, pos.z ^ 1); break;
				}

				// translation
				trans = center + search_matrix[number][dir] - face_adjacent_node_pos;

				switch (dir) {
				case 0:
					for (uint8_t y = 0; y < 2; ++y)
						for (uint8_t z = 0; z < 2; ++z) {
							neighbor[(pos.x ^ 1) + trans.x][y + trans.y][z + trans.z] = &temp->children[z << 2 | y << 1 | (pos.x ^ 1)];
						}
					break;
				case 1:
					for (uint8_t x = 0; x < 2; ++x)
						for (uint8_t z = 0; z < 2; ++z) {
							neighbor[x + trans.x][(pos.y ^ 1) + trans.y][z + trans.z] = &temp->children[z << 2 | (pos.y ^ 1) << 1 | x];
						}
					break;
				case 2:
					for (uint8_t x = 0; x < 2; ++x)
						for (uint8_t y = 0; y < 2; ++y) {
							neighbor[x + trans.x][y + trans.y][(pos.z ^ 1) + trans.z] = &temp->children[(pos.z ^ 1) << 2 | y << 1 | x];
						}
					break;
				}
			}
		}

		Position<int8_t> edge_adjacent_node_pos;
		// set 6 nodes of 27 in total
		// Find 6 cube nodes adjacent to cnode's 3 edges
		for (uint8_t dir = 0; dir < 3; ++dir) {
			temp = parent_magic_cube(center + search_matrix[number][dir + 3]);
			if (temp != nullptr) {
				if (temp->children == nullptr) temp->init_eight_children();

				switch (dir) {
				case 0:
					edge_adjacent_node_pos = Position<uint8_t>(pos.x, pos.y ^ 1, pos.z ^ 1);
					trans = center + search_matrix[number][dir + 3] - edge_adjacent_node_pos;
					neighbor[0 + trans.x][(pos.y ^ 1) + trans.y][(pos.z ^ 1) + trans.z] = &temp->children[(pos.z ^ 1) << 2 | (pos.y ^ 1) << 1 | 0];
					neighbor[1 + trans.x][(pos.y ^ 1) + trans.y][(pos.z ^ 1) + trans.z] = &temp->children[(pos.z ^ 1) << 2 | (pos.y ^ 1) << 1 | 1];
					break;
				case 1:
					edge_adjacent_node_pos = Position<uint8_t>(pos.x ^ 1, pos.y, pos.z ^ 1);
					trans = center + search_matrix[number][dir + 3] - edge_adjacent_node_pos;
					neighbor[(pos.x ^ 1) + trans.x][0 + trans.y][(pos.z ^ 1) + trans.z] = &temp->children[(pos.z ^ 1) << 2 | 0 << 1 | (pos.x ^ 1)];
					neighbor[(pos.x ^ 1) + trans.x][1 + trans.y][(pos.z ^ 1) + trans.z] = &temp->children[(pos.z ^ 1) << 2 | 1 << 1 | (pos.x ^ 1)];
					break;
				case 2:
					edge_adjacent_node_pos = Position<uint8_t>(pos.x ^ 1, pos.y ^ 1, pos.z);
					trans = center + search_matrix[number][dir + 3] - edge_adjacent_node_pos;
					neighbor[(pos.x ^ 1) + trans.x][(pos.y ^ 1) + trans.y][0 + trans.z] = &temp->children[0 << 2 | (pos.y ^ 1) << 1 | (pos.x ^ 1)];
					neighbor[(pos.x ^ 1) + trans.x][(pos.y ^ 1) + trans.y][1 + trans.z] = &temp->children[1 << 2 | (pos.y ^ 1) << 1 | (pos.x ^ 1)];
					break;
				}
			}
		}

		// set 1 node of 27 in total
		// Find 1 cube node adjacent to cnode's 1 vertex
		Position<int8_t> vertex_adjacent_node_pos(pos.x ^ 1, pos.y ^ 1, pos.z ^ 1);
		temp = parent_magic_cube(center + search_matrix[number][6]);
		if (temp != nullptr) {
			if (temp->children == nullptr) temp->init_eight_children();
			trans = center + search_matrix[number][6] - vertex_adjacent_node_pos;
			neighbor[(pos.x ^ 1) + trans.x][(pos.y ^ 1) + trans.y][(pos.z ^ 1) + trans.z] = &temp->children[(pos.z ^ 1) << 2 | (pos.y ^ 1) << 1 | (pos.x ^ 1)];
		}
	}

}
