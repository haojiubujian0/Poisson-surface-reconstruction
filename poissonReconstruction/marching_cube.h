/**
  * FileName: marching_cube.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date:
  * Description: 
***/

#ifndef  POISSON_RECONSTRUCTION_MARCHING_CUBE_H
#define POISSON_RECONSTRUCTION_MARCHING_CUBE_H
#include <iostream>
#include <vector>

namespace poisson_reconstruction {
	class Marching_cube {
	public:
		static const int edgeMask[1 << 8];
		static const int triangles[1 << 8][3 * 5 + 1];
		static uint8_t get_cube_index(const double corners[8], double isovalue);
		static int add_triangles(int cube_index, std::vector<int>& tris);
	};

}

#endif // ! POISSON_RECONSTRUCTION_MARCHING_CUBE_H
