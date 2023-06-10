/**
  * FileName: polynomial.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date: 2023/3/17/20:36
  * Description: A simple class for representing polynomial functions.
***/

#ifndef  POISSON_RECONSTRUCTION_POLYNOMIAL_H
#define  POISSON_RECONSTRUCTION_POLYNOMIAL_H

#include <cstring>
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <iostream>

namespace poisson_reconstruction {
	// A polynomial defined by element type T and degree D
	template <uint16_t D, typename T = double>
	class Polynomial {
		static const double epsilon;
	public:
		Polynomial();
		~Polynomial() = default;
		
		template <uint16_t Drhs, typename Trhs>
		Polynomial(const Polynomial<Drhs, Trhs> &rhs) {
			std::memset(coefficients, 0, sizeof(T) * (D + 1));
			for (int i = 0; i <= std::min(D, Drhs); ++i) {
				coefficients[i] = rhs[i];
			}
		}

		template <uint16_t Drhs, typename Trhs>
		Polynomial& operator=(const Polynomial<Drhs, Trhs> &rhs) {
			// handle self-assignment
			if (D == Drhs && std::is_same<T, Trhs>::value) {
				if ((void*)&rhs == (void*)this) return *this;
			}
			
			memset(coefficients, 0, sizeof(T) * (D + 1));
			for (int i = 0; i <= std::min(D, Drhs); ++i) {
				coefficients[i] = rhs[i];
			}
			return *this;
		}

		template <uint16_t Drhs, typename Trhs>
		Polynomial& operator+=(const Polynomial<Drhs, Trhs>& rhs) {
			for (int i = 0; i <= std::min(D, Drhs); ++i) {
				coefficients[i] += T(rhs[i]);
			}
			return *this;
		}

		uint16_t degree() const { return D; }

		// Get the coefficient of degree i, i.e x^i
		T& operator[](unsigned short i) { return coefficients[i]; }
		const T& operator[](unsigned short i) const { return coefficients[i]; }
		
		// Translation operation, i.e P(x+scalar)
		Polynomial translation(T scalar) const ;
		
		// Scaling operation, i.e P(x/scalar)
		Polynomial scale(T scalar) const ;

		// P(x) = c0*x^0 + c1*x^1 + ... + ci*x^i + ... + cD*x^D
		// ans  = ... + ci/(i+1)*x^(i+1) + ...    
		Polynomial<D + 1, T> antiderivative() const {
			Polynomial<D + 1, T> ans;
			for (int i = 0; i <= D; ++i) {
				// handle the i-th formular
				ans[i + 1] = coefficients[i] / (i + 1);
			}
			return ans;
		}

		// derived function
		// assume that D >= 1
		// P(x) = c0*x^0 + c1*x^1 + ... + ci*x^i + ... + cD*x^D
		// ans  = ... + ci*i*x^(i-1) + ...   
		Polynomial<D - 1, T> derivative() const {
			Polynomial<D - 1, T> ans;
			// the 0-degree is discarded
			for (int i = 1; i <= D; ++i) {
				ans[i - 1] = coefficients[i] * i;
			}
			return ans;
		}

		// the definite integral from s to e
		T integral(T s, T e) const {
			if (std::fabs(e - s) < epsilon) return 0;
			// primitive function 
			auto primitive = antiderivative();
			return primitive(e) - primitive(s);
		}

		// count the function value, i.e P(s)
		// P(s) = c0*(s)^0 + c1*(s)^1 + ... + ci*(s)^i + ... + cD*(s)^D 
		T operator()(T s) const {
			T ans = 0, temp = 1;
			for (int i = 0; i <= D; ++i) {
				ans += coefficients[i] * temp;
				temp *= s;
			}
			return ans;
		}

	private:
		T coefficients[D + 1];
	};

	template<uint16_t D, typename T>
	const double Polynomial<D, T>::epsilon = 1e-15;

	template <uint16_t D, typename T>
	Polynomial<D, T>::Polynomial() {
		std::memset(coefficients, 0, sizeof(T) * (D + 1));
	}

	// P(x + scalar) = c0*(x+scalar)^0 + c1*(x+scalar)^1 + ... + ci*(x+scalar)^i + ... + cD*(x+scalar)^D 
	template<uint16_t D, typename T>
	Polynomial<D, T> Polynomial<D, T>::translation(T scalar) const {
		Polynomial ans;
		T temp = 1;
		double C;
		for (int i = 0; i <= D; ++i) {
			// handle the i-th expanded formular, i.e (x+scalar)^i 
			// C{i, k} = i! / ((i-k)!k!)  C{i, k-1} = i! / ((i-k+1)!(k-1)!)
			// C{i, k} = C{i, k-1} * (i-k+1) / k
			temp = 1;
			C = 1;
			for (int k = 0; k <= i; ++k) {
				ans[i-k] += coefficients[i] * temp * C;
				temp *= scalar;
				C *= (double)(i - (k+1) + 1) / (k+1);
			}
		}
		return ans;
	}

	// P(x/scalar) = c0*(x/scalar)^0 + c1*(x/scalar)^1 + ... + ci*(x/scalar)^i + ... + cD*(x/scalar)^D
	template<uint16_t D, typename T>
	Polynomial<D, T> Polynomial<D, T>::scale(T scalar) const {
		Polynomial ans;
		T temp = 1;
		for (int i = 0; i <= D; ++i) {
			// handle the i-th expanded formular, i.e (x/scalar)^i
			ans[i] = coefficients[i] * (1.0/temp);
			temp *= scalar;
		}
		return ans;
	}

	//---------------------------------------------------------------------------------------------------
	// non-member functions
	//---------------------------------------------------------------------------------------------------

	// Polynomial multiplication
	template <uint16_t Dlhs, typename Tlhs,
			  uint16_t Drhs, typename Trhs>
	auto operator*(const Polynomial<Dlhs, Tlhs> &lhs, const Polynomial<Drhs, Trhs> &rhs) {
		Polynomial <Dlhs + Drhs, decltype(Tlhs() + Trhs())> ans;
		for(int i = 0; i <= Dlhs; ++i)
			for (int j = 0; j <= Drhs; ++j) {
				ans[i + j] += lhs[i] * rhs[j];
			}
		return ans;
	}

	// Polynomial*scalar multiplication 
	template <uint16_t D, typename T>
	Polynomial<D, T> operator*(const Polynomial<D, T>& lhs, T scalar) {
		Polynomial<D, T> ans(lhs);
		for (int i = 0; i <= D; ++i) {
			ans[i] *= scalar;
		}
		return ans;
	}

	// Polynomial*scalar division 
	template <uint16_t D, typename T, typename T2>
	Polynomial<D, T> operator/(const Polynomial<D, T>& lhs, T2 scalar) {
		Polynomial<D, T> ans(lhs);
		for (int i = 0; i <= D; ++i) {
			ans[i] /= (T)scalar;
		}
		return ans;
	}

	// Polynomial addition
	template <uint16_t Dlhs, typename Tlhs,
			  uint16_t Drhs, typename Trhs>
	auto operator+(const Polynomial<Dlhs, Tlhs>& lhs, const Polynomial<Drhs, Trhs>& rhs) {
		unsigned short D = std::max(Dlhs, Drhs);
		decltype(Tlhs() + Trhs()) temp = 0;
		Polynomial <D, decltype(temp)> ans;
		for (int i = 0; i <= std::max(Dlhs, Drhs); ++i) {
			temp = 0;
			if (i <= Dlhs) temp += lhs[i];
			if (i <= Drhs) temp += rhs[i];
			ans[i] = temp;
		}
		return ans;
	}

	// Polynomial+scalar addition 
	template <uint16_t D, typename T>
	Polynomial<D, T> operator+(const Polynomial<D, T>& lhs, T scalar) {
		Polynomial<D, T> ans(lhs);
		ans[0] += scalar;
		return ans;
	}
	
	// Polynomial subtraction
	template <uint16_t Dlhs, typename Tlhs,
			  uint16_t Drhs, typename Trhs>
	auto operator-(const Polynomial<Dlhs, Tlhs>& lhs, const Polynomial<Drhs, Trhs>& rhs) {
		unsigned short D = std::max(Dlhs, Drhs);
		decltype(Tlhs() - Trhs()) temp = 0;
		Polynomial <D, decltype(temp)> ans;
		for (int i = 0; i <= std::max(Dlhs, Drhs); ++i) {
			temp = 0;
			if (i <= Dlhs) temp = lhs[i];
			if (i <= Drhs) temp -= rhs[i];
			ans[i] = temp;
		}
		return ans;
	}

	// Polynomial-scalar subtraction
	template <uint16_t D, typename T>
	Polynomial<D, T> operator-(const Polynomial<D, T> &lhs, T scalar) {
		Polynomial<D, T> ans(lhs);
		ans[0] -= scalar;
		return ans;
	}

	// Polynomial-scalar subtraction
	template <uint16_t D, typename T>
	Polynomial<D, T> operator-(T scalar, const Polynomial<D, T>& lhs) {
		Polynomial<D, T> ans(lhs);
		ans[0] = scalar - ans[0];
		for (int i = 1; i <= D; ++i) {
			ans[i] = -ans[i];
		}
		return ans;
	}

	template<uint16_t D, typename T>
	std::ostream& operator<<(std::ostream& os, const Polynomial<D, T>& poly) {
		for (int i = 0; i <= D; ++i) {
			os << poly[i] << "*x^" << i;
			if (i < D && poly[i+1] >= 0) os << "+";
		}
		return os;
	}
}
#endif // ! POISSON_RECONSTRUCTION_POLYNOMIAL_H
