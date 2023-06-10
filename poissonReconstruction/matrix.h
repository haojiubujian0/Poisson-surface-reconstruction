/**
  * FileName: matrix.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date: 
  * Description: Implement sparsely symmetric matrix
***/

#ifndef POISSON_RECONSTRUCTION_MATRIX_H
#define POISSON_RECONSTRUCTION_MATRIX_H

#include <iostream>
#include <vector>
#include <utility>

namespace poisson_reconstruction {
	template<typename T>
	class Entry {
	public:
		uint32_t col{ 0 };
		T value{ 0 };
	};

	template<typename T>
	class Sparse_sym_matrix {
	public:
		Sparse_sym_matrix(uint32_t rows):matrix(rows) { }
		Sparse_sym_matrix() = default;

		uint32_t get_rows() const { return matrix.size(); }
		uint32_t get_row_size(uint32_t r) const { return matrix[r].size(); }

		void push_back(uint32_t r, const Entry<T>& entry);
		void push_back(std::vector<Entry<T>>&& row_elemens) {
			matrix.push_back(std::move(row_elemens));
		}

		std::vector<Entry<T>>& operator[](uint32_t r) {
			return matrix[r];
		}
		const std::vector<Entry<T>>& operator[](uint32_t r) const {
			return matrix[r];
		}

	private:
		std::vector<std::vector<Entry<T>>> matrix;
	};

	template<typename T>
	inline void Sparse_sym_matrix<T>::push_back(uint32_t r, const Entry<T>& entry) {
		matrix[r].push_back(entry);
	}

	//---------------------------------------------------------------------------------------------------
	// non-member functions
	//---------------------------------------------------------------------------------------------------

	// assume that v1 and v2 have same size
	// The dot product of the vectors
	template<typename T1, typename T2>
	decltype(T1() * T2()) operator*(const std::vector<T1>& v1, const std::vector<T2>& v2) {
		decltype(T1()* T2()) ans = 0;

		uint32_t cols = v1.size();
		for (uint32_t i = 0; i < cols; ++i) {
			ans += v1[i] * v2[i];
		}
		
		return ans;
	}

	// assume that v1 and v2 have same size
	// Vector subtraction
	template<typename T1, typename T2>
	std::vector<decltype(T1() - T2())> operator-(const std::vector<T1>& v1, const std::vector<T2>& v2) {
		std::vector<decltype(T1() - T2())> ans(v1.size());

		uint32_t cols = v1.size();
		for (uint32_t i = 0; i < cols; ++i) {
			ans[i] = v1[i] - v2[i];
		}

		return ans;
	}

	// assume that v1 and v2 have same size
	// vector addition
	template<typename T1, typename T2>
	std::vector<decltype(T1() + T2())> operator+(const std::vector<T1>& v1, const std::vector<T2>& v2) {
		std::vector<decltype(T1() + T2())> ans(v1.size());

		uint32_t cols = v1.size();
		for (uint32_t i = 0; i < cols; ++i) {
			ans[i] = v1[i] + v2[i];
		}

		return ans;
	}

	// vector-scalar multiplication
	template<typename T1, typename T2>
	std::vector<decltype(T1()* T2())> operator*(const std::vector<T1>& v, T2 s) {
		std::vector<decltype(T1()* T2())> ans(v.size(), 0);

		uint32_t cols = v.size();
		for (uint32_t i = 0; i < cols; ++i) {
			ans[i] = v[i] * s;
		}

		return ans;
	}

	// assume that matrix: m*n, v: n*1
	template<typename T1, typename T2>
	std::vector<decltype(T1()* T2())> operator*(const Sparse_sym_matrix<T1>& matrix, const std::vector<T2>& v) {
		uint32_t rows = matrix.get_rows();
		std::vector<decltype(T1() * T2())> ans(rows, 0);

		uint32_t cols;
		uint32_t real_col_num;
		for (uint32_t i = 0; i < rows; ++i) {
			cols = matrix.get_row_size(i);
			for (uint32_t j = 0; j < cols; ++j) {
				real_col_num = matrix[i][j].col;

				ans[i] += matrix[i][j].value * v[real_col_num];
			}
		}

		return ans;
	}

	// conjugate gradient algorithm
	// matrix * x = b, parameter solution is an approximate solution 
	// matrix: m*m
	template<typename T1, typename T2>
	void CG(const Sparse_sym_matrix<T1>& matrix, const std::vector<T2>& b, std::vector<decltype(T1()*T2())>& solution, 
			uint32_t max_iters, double epsilon = 1e-6) 
	{
		std::cout << "max iters: " << max_iters << std::endl;
		solution.clear();
		solution.resize(b.size(), 0);

		double alpha = 0, beta = 0;
		std::vector<decltype(T1()* T2())> r(b.size(), 0);
		std::vector<decltype(T1()* T2())> p(b.size(), 0);

		// set r0 and p0
		r = b - matrix * solution;
		p = r;

		std::vector<decltype(T1()* T2())> temp(b.size(), 0);
		std::vector<decltype(T1()* T2())> r_next(b.size(), 0);
		double rr = 0, tp = 0;
		uint32_t k;
		for (k = 0; k < max_iters; ++k) {
			std::cout << "curr iter: " << k << std::endl;
			temp = matrix * p;
			rr = r * r;
			tp = temp * p;
			if (std::fabs(rr) < epsilon || std::fabs(tp) < epsilon) break;
			alpha = rr / tp;
			solution = solution + p * alpha;
			r_next = r - temp * alpha;
			beta = (r_next * r_next) / (r * r);
			r.swap(r_next);
			p = r + p * beta;
		}
		std::cout << "iteration: " << k << std::endl;
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& os, const Sparse_sym_matrix<T>& matrix) {
		int rows = matrix.get_rows();
		os << "sparsely symmetric matrix: (real col, value), rows " << rows << '\n';
		int cols;
		for (int i = 0; i < rows; ++i) {
			cols = matrix.get_row_size(i);
			os << "row " << i << ": \n";
			for (int j = 0; j < cols; ++j) {
				os << "( " << matrix[i][j].col << " , " << matrix[i][j].value << ")";
			}
			os << '\n';
		}
		return os;
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
		int num = v.size();
		os << "vector: rows " << num << '\n';
		for (int i = 0; i < num; ++i) {
			os << v[i] << " ";
		}
		return os;
	}
}

#endif // !POISSON_RECONSTRUCTION_MATRIX_H