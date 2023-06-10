/**
  * FileName: gaussian.h
  * Author: Zeng Zhenxiang
  * Version: 1.0
  * Date: 2023/3/18/16:22
  * Description: Simulating Gaussian filters using convolution of box filters.
    i.e: 
        Box(t) = 1/(2*radius)  |t| < radius
               = 0             otherwise
        The default setting: radius = 0.5, so that
        Box(t) = 1  |t| < 0.5
               = 0  otherwise
        Gaussian(t) = Box*n , the n-th convolution of a box filter with itself   
***/

#ifndef POISSON_RECONSTRUCTION_GAUSSIAN_H
#define POISSON_RECONSTRUCTION_GAUSSIAN_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include "polynomial.h"

namespace poisson_reconstruction {
    
    // Polynomial element, i.e a polynomial with a start
    // A Polynomial element is defined in [start, inf)
    template <uint16_t D, typename T = double>
    class Poly_element {
    public:
        T start;
        Polynomial<D, T> poly;
        
        template <uint16_t Drhs, typename Trhs>
        Poly_element& operator+=(const Poly_element<Drhs, Trhs>& rhs) {
            start = std::max(start, (T)rhs.start);
            poly += rhs.poly;
            return (*this);
        }

        template <uint16_t Drhs, typename Trhs>
        Poly_element& operator=(const Poly_element<Drhs, Trhs>& rhs) {
            start = std::max(start, (T)rhs.start);
            poly = rhs.poly;
            return (*this);
        }

        // indefinite integral
        Poly_element<D + 1, T> antiderivative() const {
            Poly_element<D+1, T> ans;
            ans.start = start;
            ans.poly = poly.antiderivative();
            return ans;
        }

        // derived function
        Poly_element<D - 1, T> derivative() const {
            Poly_element<D - 1, T> ans;
            ans.start = start;
            ans.poly = poly.derivative();
            return ans;
        }

        // Function translation operation
        Poly_element translation(T s) const {
            Poly_element ans;
            ans.start = start - s;
            ans.poly = poly.translation(s);
            return ans;
        }

        // scaling operation
        // assume that s > 0
        Poly_element scale(T s) const {
            Poly_element ans;
            ans.start = start*s;
            ans.poly = poly.scale(s);
            return ans;
        }

        // count the function value, i.e P(s)
        T operator()(T s) const {
            if (s < start) return 0;
            return poly(s);
        }

        // the definite integral from s to e
        T integral(T s, T e) const {
            s = std::max(s, start);
            if (s >= e) return 0;
            return poly.integral(s, e);
        }
    };

    // A Gaussian approximation is a symmetric polynomial, 
    // but a symmetric polynomial is not always a gaussian approximation.
    template <uint16_t D, typename T = double>
    class Symmetric_poly {
        static const double epsilon;
        using p_element = Poly_element<D, T>;

        void sort_polys() {
            std::sort(polys.begin(), polys.end(), [](const p_element& lhs, const p_element& rhs) {
                return lhs.start < rhs.start;
                });
        }

        // Assume that the elements in the container polys has sorted.
        void merge() {
            int n = polys.size();
            if (n <= 1) return;
            std::vector<p_element> ans;
            p_element temp = polys[0];
            for (int i = 1; i < n; ++i) {
                if (std::fabs(polys[i].start - polys[i - 1].start) < epsilon) {
                    temp += polys[i];
                }
                else {
                    ans.push_back(temp);
                    temp = polys[i];
                }
            }
            ans.push_back(temp);
            polys.swap(ans);
        }

    public:
        Symmetric_poly() = default;
        ~Symmetric_poly() = default;

        template <uint16_t Drhs, typename Trhs>
        Symmetric_poly(const Symmetric_poly<Drhs, Trhs>& rhs) { 
            size_t num = rhs.polys.size();
            polys.resize(num);
            for (size_t i = 0; i < num; ++i) {
                polys[i] = rhs.polys[i];
            }
        }

        template <uint16_t Drhs, typename Trhs>
        Symmetric_poly<D, T>& operator=(const Symmetric_poly<Drhs, Trhs>& rhs);

        // access the i-th polynomial element.
        p_element& operator[](size_t i) { return polys[i]; }
        const p_element& operator[](size_t i) const { return polys[i]; }
        
        // add a polynomial element into the container, polys.
        void push_back(const p_element& pe) { polys.push_back(pe); sort_polys(); merge(); }
        
        size_t size() const { return polys.size(); }

        // count the function value, i.e P(s)
        T operator()(T s) const {
            T ans = 0;
            for (int i = 0; i < polys.size(); ++i) {
                ans += polys[i](s);
            }
            return ans;
        }

        // indefinite integral
        Symmetric_poly<D + 1, T> antiderivative() const {
            Symmetric_poly<D + 1, T> ans;
            for (int i = 0; i < polys.size(); ++i) {
                ans.push_back(polys[i].antiderivative());
            }
        }

        // derived function
        Symmetric_poly<D - 1, T> derivative() const {
            Symmetric_poly<D - 1, T> ans;
            for (int i = 0; i < polys.size(); ++i) {
                ans.push_back(polys[i].derivative());
            }
            return ans;
        }

        // Function translation operation
        Symmetric_poly translation(T s) const {
            Symmetric_poly ans;
            for (int i = 0; i < polys.size(); ++i) {
                ans.push_back(polys[i].translation(s));
            }
            return ans;
        }

        // scaling operation 
        // assume that s > 0
        Symmetric_poly scale(T s) const {
            Symmetric_poly ans;
            for (int i = 0; i < polys.size(); ++i) {
                ans.push_back(polys[i].scale(s));
            }
            return ans;
        }

        // the definite integral from s to e
        T integral(T s, T e) const {
            T ans = 0;
            for (int i = 0; i < polys.size(); ++i) {
                ans += polys[i].integral(s, e);
            }
            return ans;
        }

        const std::vector<p_element>& get_polys() const { return polys; }
        
    //private:
        std::vector<p_element> polys;
    };
    
    template<uint16_t D, typename T>
    const double Symmetric_poly<D, T>::epsilon = 1e-15;

    template<uint16_t D, typename T>
    template<uint16_t Drhs, typename Trhs>
    inline Symmetric_poly<D, T>& Symmetric_poly<D, T>::operator=(const Symmetric_poly<Drhs, Trhs>& rhs) {
        // handle self-assignment
        if (this == &rhs) return *this;

        polys = rhs.polys;
        return *this;
    }

    //---------------------------------------------------------------------------------------------------
    // non-member functions
    //---------------------------------------------------------------------------------------------------
    
    // Polynomial element multiplication
    template <uint16_t Dlhs, typename Tlhs,
              uint16_t Drhs, typename Trhs>
    inline auto operator*(const Poly_element<Dlhs, Tlhs>& lhs, const Poly_element<Drhs, Trhs>& rhs) {
        Poly_element<Dlhs + Drhs, decltype(Tlhs() + Trhs())> ans;
        ans.start = std::max(lhs.start, rhs.start);
        ans.poly = lhs.poly * rhs.poly;
        return ans;
    }

    // Symmetric Polynomial multiplication
    template <uint16_t Dlhs, typename Tlhs,
              uint16_t Drhs, typename Trhs>
    inline auto operator*(const Symmetric_poly<Dlhs, Tlhs>& lhs, const Symmetric_poly<Drhs, Trhs>& rhs) {
        Symmetric_poly<Dlhs + Drhs, decltype(Tlhs() + Trhs())> ans;

        for(int i = 0; i < lhs.size(); ++i)
            for (int j = 0; j < rhs.size(); ++j) {
                // handle lhs[i] * rhs[j]
                ans.push_back(lhs[i] * rhs[j]);
            }

        return ans;
    }

    // Symmetric Polynomial division
    template <uint16_t D, typename T, typename T2>
    inline auto operator/(const Symmetric_poly<D, T>& lhs, T2 scalar) {
        Symmetric_poly<D, T> ans;
        Poly_element<D, T> temp;

        int num = lhs.size();
        for (int i = 0; i < num; ++i) {
            temp.start = lhs[i].start;
            temp.poly = lhs[i].poly / scalar;
            ans.push_back(temp);
        }

        return ans;
    }

    // A gaussian approximation representing the N-th convolution of a box filter with itself.
    template <uint16_t N, typename T>
    Symmetric_poly<N-1, T> Gaussian_approximation(T radius = 0.5) {
        auto sp = Gaussian_approximation<N - 1, T>(radius);
        // when call this function, N is not equal to 1
        
        Symmetric_poly<N - 1, T> ans;
        Poly_element<N - 1, T> temp;
        Polynomial<N - 1, T> primitive;

        // the convolution between Box*(n-1) and Box
        for (int i = 0; i < sp.size(); ++i) {
            // handle the i-th polynomial element
            primitive = sp[i].poly.antiderivative();

            temp.start = sp[i].start - radius;
            temp.poly = primitive.translation(radius) - primitive(sp[i].start);
            ans.push_back(temp);
            
            temp.start = sp[i].start + radius;
            temp.poly = primitive(sp[i].start) - primitive.translation(-radius);
            ans.push_back(temp);
        }
        return ans;
    }

    // box filter
    // B(t) = 1/(2*radius) |t| < radius
    //      = 0           otherwise
    template <>
    Symmetric_poly<(uint16_t)0, double> Gaussian_approximation<(uint16_t)1, double>(double radius) {
        Symmetric_poly<0, double> box;
        Poly_element<0, double> p;
        p.start = -1 * radius; p.poly[0] = 1 / (2 * radius);
        box.push_back(p);
        p.start = radius; p.poly[0] = -1 / (2 * radius);
        box.push_back(p);
        return box;
    }

    template <>
    Symmetric_poly<(uint16_t)0, float> Gaussian_approximation<(uint16_t)1, float>(float radius) {
        Symmetric_poly<0, float> box;
        Poly_element<0, float> p;
        p.start = -1 * radius; p.poly[0] = 1 / (2 * radius);
        box.push_back(p);
        p.start = radius; p.poly[0] = -1 / (2 * radius);
        box.push_back(p);
        return box;
    }

    template<uint16_t D, typename T>
    std::ostream& operator<<(std::ostream& os, const Symmetric_poly<D, T>& sym_poly) {
        Polynomial <D, T> poly_temp;
        int num = sym_poly.size();
        for (int i = 0; i < num; ++i) {
            os << "[ " << sym_poly[i].start << " , ";
            if (i + 1 < num) os << sym_poly[i + 1].start << " ]";
            else os << "infinity )";
            os << '\t';
            poly_temp += sym_poly[i].poly;
            os << poly_temp << std::endl;
        }
        return os;
    }
}

#endif