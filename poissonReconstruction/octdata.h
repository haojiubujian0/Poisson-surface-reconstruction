/**
  * FileName: octdata.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date: 2023/3/19/11:00
  * Description: the function data used in octree.
    Some definitions are as follows: for a cube node o, o.center = (o.cx, o.cy, o.cz), o.width
		Bo(x) = Box*N((x-o.cx)/o.width)	(convolution)
		Bo(y) = Box*N((y-o.cy)/o.width)	(convolution)
		Bo(z) = Box*N((z-o.cz)/o.width)	(convolution)
		Fo(x, y, z) = Bo(x)*Bo(y)*Bo(z)*(1/o.width^3)
***/

#ifndef POISSON_RECONSTRUCTION_OCTDATA_H
#define POISSON_RECONSTRUCTION_OCTDATA_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "gaussian.h"

namespace poisson_reconstruction {
	
	//---------------------------------------------------------------------------------------------------
	// class Function_set
	//---------------------------------------------------------------------------------------------------
	
	// Stores all the functions with all different resolutions and offsets for max depth.
	// The x, y and z directions are the same
	template <uint16_t N = 3, typename T = double>
	class Function_set {
	public:
		Function_set() = default;
		void set_function_set(uint8_t max_d);
		
		static void get_depth_offset(uint16_t number, uint8_t& depth, uint16_t& offset);

		const Symmetric_poly<N - 1, T>& get_base_func() const { return base_func; }
		uint16_t get_total_num() const { return total_num; }

		// number: the cube node number in x, y or z direction
		// this cube number is transformed from the cube node's Offset
		// Offset.x  -> the number of the node in the x direction, y z is in a similar way
		// Return the function Bo(x), Bo(y) or Bo(z)
		Symmetric_poly<N - 1, T>& operator[](uint16_t number) { return functions[number]; }
		const Symmetric_poly<N - 1, T>& operator[](uint16_t number) const { return functions[number]; }

	private:
		// The max depth of octree
		uint8_t max_depth{ 0 };
		uint16_t total_num{ 0 };
		// Box*N, i.e the n-th convolution of a box filter with itself  
		Symmetric_poly<N - 1, T> base_func;
		// All the functions that are scaled and shifted by the basis function
		std::vector<Symmetric_poly<N - 1, T>> functions;
	};

	template<uint16_t N, typename T>
	inline void Function_set<N, T>::set_function_set(uint8_t max_d) {

		//std::ofstream polys_file;
		//polys_file.open("polys.txt", polys_file.out | polys_file.trunc);

		//if (!polys_file) {
		//	throw std::runtime_error("can't open polys.txt");
		//}

		max_depth = max_d;
		if (max_depth < 0) return;
		total_num = (1 << (max_depth + 1)) - 1;
		base_func = Gaussian_approximation<N, T>(0.5);
		base_func = base_func / base_func(0);

		//polys_file << "baseFunction: \n" << base_func;
		//polys_file << "dBaseFunction: \n" << base_func.derivative();

		functions.reserve(total_num);
		uint8_t depth;
		uint16_t offset;
		double width, center;
		for (uint16_t number = 0; number < total_num; ++number) {
			get_depth_offset(number, depth, offset);
			width = 1.0 / (1 << depth);
			center = (offset + 0.5) * width;
			functions.push_back(base_func.scale(width).translation(-center));

			//polys_file << "µÚ " << (int)number << " ¸ö" << std::endl;
			//polys_file << "w = " << width << " , c = " << center << " genereate: \n" << functions[number];
		}
		//polys_file.close();
	}

	// here: number = offset + (1<<depth) - 1
	// offset and depth are properties of octree node
	// offset can be Offset.x , Offset.y or Offset.z
	template<uint16_t N, typename T>
	inline void Function_set<N, T>::get_depth_offset(uint16_t number, uint8_t& depth, uint16_t& offset) {
		int i = number + 1;
		depth = 0;
		while (i != 1) {
			depth += 1;
			i /= 2;
		}
		// inverse transformation
		offset = number + 1 - (1 << depth);
	}

	//---------------------------------------------------------------------------------------------------
	// class Octdata
	//---------------------------------------------------------------------------------------------------

	template <uint16_t N = 3, typename T = double>
	class Octdata {
	public:
		Octdata(uint8_t max_d): max_depth(max_d), sampling_interval(1.0/(1<<(max_d+1))) { func_set.set_function_set(max_d); }

		void set_ff_dot();
		void set_df_dot();
		void set_dd_dot();
		void set_values();

		const std::vector<std::vector<T>>& get_ff_dot() { return ff_dot; }
		const std::vector<std::vector<T>>& get_df_dot() { return df_dot; }
		const std::vector<std::vector<T>>& get_dd_dot() { return dd_dot; }
		const std::vector<std::vector<T>>& get_values() { return values; }

		const Function_set<N, T>& get_func_set() const { return func_set; }

		void write_to_disk(const std::string& filename, const std::vector<std::vector<T>>& table);

	private:
		uint8_t max_depth;
		// the radius of Box*N
		T radius{N*0.5};

		Function_set<N, T> func_set;
		
		// dot product: <Bo, Bo2>
		// this is: < func_set.functions[i] , func_set.functions[j] > 
		std::vector<std::vector<T>> ff_dot;
		
		// dot product: <dBo, Bo2>
		// this is: < d func_set.functions[i] , func_set.functions[j] >
		std::vector<std::vector<T>> df_dot;

		// dot product: <d(dBo), Bo2>
		// this is: < dd func_set.functions[i] , func_set.functions[j] >
		std::vector<std::vector<T>> dd_dot;

		// functions values
		// values[i][j] means the value of func_set.functions[i] in the coordinate j/(1<<max_depth)
		T sampling_interval;  // 1/(1<<max_depth)
		std::vector<std::vector<T>> values;
	};

	template<uint16_t N, typename T>
	inline void Octdata<N, T>::set_ff_dot() {
		if (!ff_dot.empty()) return;

		uint16_t nums = func_set.get_total_num();
		ff_dot = std::vector<std::vector<T>>(nums, std::vector<T>(nums));

		const auto& base_func = func_set.get_base_func();
		T temp = 0;
		uint16_t off1, off2;
		uint8_t d1, d2;
		T w1, c1, w2, c2;
		for(int i = 0; i < nums; ++i)
			for (int j = 0; j <= i; ++j) {
				// handle the dot product: < func_set.functions[i] , func_set.functions[j] >

				// get the func: Boi(t) = Box*N((t-c1)/w1)
				Function_set<N, T>::get_depth_offset(i, d1, off1);
				w1 = 1.0 / (1 << d1); 
				c1 = (off1 + 0.5) * w1;

				// get the func: Boj(t) = Box*N((t-c2)/w2)
				Function_set<N, T>::get_depth_offset(j, d2, off2);
				w2 = 1.0 / (1 << d2);
				c2 = (off2 + 0.5) * w2;

				// <Boi, Boj>
				temp = (base_func * base_func.scale(w2 / w1).translation(-((c2 - c1) / w1))).integral(-radius, radius) * w1;
				
				// symmetry
				ff_dot[i][j] = ff_dot[j][i] = temp;
			}
	}

	template<uint16_t N, typename T>
	inline void Octdata<N, T>::set_df_dot() {
		if (!df_dot.empty()) return;

		uint16_t nums = func_set.get_total_num();
		df_dot = std::vector<std::vector<T>>(nums, std::vector<T>(nums));

		const auto& base_func = func_set.get_base_func();
		T temp = 0;
		uint16_t off1, off2;
		uint8_t d1, d2;
		T w1, c1, w2, c2;
		for (int i = 0; i < nums; ++i)
			for (int j = 0; j <= i; ++j) {
				// handle the dot product: < d func_set.functions[i] , func_set.functions[j] >  

				// get the func: Boi(t) = Box*N((t-c1)/w1)
				Function_set<N, T>::get_depth_offset(i, d1, off1);
				w1 = 1.0 / (1 << d1);
				c1 = (off1 + 0.5) * w1;

				// get the func: Boj(t) = Box*N((t-c2)/w2)
				Function_set<N, T>::get_depth_offset(j, d2, off2);
				w2 = 1.0 / (1 << d2);
				c2 = (off2 + 0.5) * w2;

				// <d Boi, Boj>
				temp = (base_func.derivative() * base_func.scale(w2 / w1).translation(-((c2 - c1) / w1))).integral(-radius, radius);

				// <d Boi, Boj> = - <d Boj, Boi>
				df_dot[i][j] = temp;
				//std::cout << temp << std::endl;
				if (i != j) df_dot[j][i] = -temp;
			}
	}

	template<uint16_t N, typename T>
	inline void Octdata<N, T>::set_dd_dot() {
		if (!dd_dot.empty()) return;

		uint16_t nums = func_set.get_total_num();
		dd_dot = std::vector<std::vector<T>>(nums, std::vector<T>(nums));

		const auto& base_func = func_set.get_base_func();

		T temp = 0;
		uint16_t off1, off2;
		uint8_t d1, d2;
		T w1, c1, w2, c2;
		for (int i = 0; i < nums; ++i)
			for (int j = 0; j <= i; ++j) {
				// handle the dot product: < dd func_set.functions[i] , func_set.functions[j] >  

				// get the func: Boi(t) = Box*N((t-c1)/w1)
				Function_set<N, T>::get_depth_offset(i, d1, off1);
				w1 = 1.0 / (1 << d1);
				c1 = (off1 + 0.5) * w1;

				// get the func: Boj(t) = Box*N((t-c2)/w2)
				Function_set<N, T>::get_depth_offset(j, d2, off2);
				w2 = 1.0 / (1 << d2);
				c2 = (off2 + 0.5) * w2;

				// <dd Boi, Boj>
				temp = (base_func.derivative().derivative() * 
						base_func.scale(w2 / w1).translation(-((c2 - c1) / w1))).integral(-radius, radius) / w1;

				// symmetry: <dd Boi, Boj> = <dd Boj, Boi>
				dd_dot[i][j] = dd_dot[j][i] = temp;
			}
	}

	template<uint16_t N, typename T>
	inline void Octdata<N, T>::set_values() {
		if (!values.empty()) return;
		
		uint16_t rows = func_set.get_total_num();
		uint16_t cols = (1 << (max_depth+1)) + 1;
		values = std::vector<std::vector<T>>(rows, std::vector<T>(cols));

		T x = 0;
		for (int i = 0; i < rows; ++i) {
			const auto& func = func_set[i];
			for (int j = 0; j < cols; ++j) {
				x = j * sampling_interval;

				values[i][j] = func(x);
				//std::cout << values[i][j] << " ";
			}
			//std::cout << std::endl;
		}
	}

	template<uint16_t N, typename T>
	inline void Octdata<N, T>::write_to_disk(const std::string& filename, const std::vector<std::vector<T>>& table) {
		std::ofstream file;
		file.open(filename, file.out | file.trunc);

		if (table.size() == 0) return;

		int rows = table.size();
		int cols = table[0].size();
		file << "rows: " << table.size() << " , " << "cols: " << table[0].size() << std::endl;
		for (int i = 0; i < rows; ++i) {
			file << "row " << i << ": ";
			for (int j = 0; j < cols; ++j) {
				file << table[i][j] << " ";
			}
			file << std::endl;
		}
		file.close();
	}

}
#endif // !POISSON_RECONSTRUCTION_OCTDATA_H
