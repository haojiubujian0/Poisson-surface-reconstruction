/**
  * FileName: mesh.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date:
  * Description:
***/

#ifndef POISSON_RECONSTRUCTION_MESH_H
#define POISSON_RECONSTRUCTION_MESH_H

#include <vector>

namespace poisson_reconstruction {
	//---------------------------------------------------------------------------------------
	// class Position
	//---------------------------------------------------------------------------------------

	template <typename T>
	struct Position {
		T x, y, z;
		Position(T xx = 0, T yy = 0, T zz = 0) : x(xx), y(yy), z(zz) { }
		template <typename Trhs>
		Position& operator=(const Position<Trhs>& rhs);
	};

	template<typename T>
	template<typename Trhs>
	inline Position<T>& Position<T>::operator=(const Position<Trhs>& rhs) {
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		return *this;
	}

	template <typename T1, typename T2>
	Position<decltype(T1() + T2())> operator+(const Position<T1>& lhs, const Position<T2>& rhs) {
		return Position<decltype(T1() + T2())>(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
	}

	template <typename T1, typename T2>
	Position<decltype(T1() - T2())> operator-(const Position<T1>& lhs, const Position<T2>& rhs) {
		return Position<decltype(T1() - T2())>(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
	}

	//---------------------------------------------------------------------------------------
	// class Mesh_data
	//---------------------------------------------------------------------------------------

	class Mesh_data {
	public:
		std::vector<Position<float>> intersection_points;  // save all the intersections
		std::vector<std::vector<int>> triangles;  // each line has 3 int
	};
}

#endif // !POISSON_RECONSTRUCTION_MESH_H
