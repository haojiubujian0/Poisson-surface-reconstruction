/**
  * FileName: magic_cube.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date:
  * Description:
***/

#ifndef POISSON_RECONSTRUCTION_MAGIC_CUBE_H
#define POISSON_RECONSTRUCTION_MAGIC_CUBE_H

#include "node.h"
#include "mesh.h"

namespace poisson_reconstruction {

	// A Rubik's cube structure
	class Magic_cube {
	public:
		Magic_cube();

		template<typename TM>
		Node* operator()(const Position<TM>& pos) { return neighbor[pos.x][pos.y][pos.z]; }

		Node* get_center_node() { return neighbor[1][1][1]; }
		// Construct the Rubik's cube structure with cnode as the center
		void construct(Node* cnode);
	private:
		Node* neighbor[3][3][3];  // means neighbor[x][y][z]
	};

}

#endif // !POISSON_RECONSTRUCTION_MAGIC_CUBE_H
