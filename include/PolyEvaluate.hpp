//
// Created by chris on 9/4/18.
//

#ifndef PROYECTO1_POLYEVALUATE_HPP
#define PROYECTO1_POLYEVALUATE_HPP

#include <iostream>
#include <string>

#include <boost/math/tools/polynomial.hpp>
#include <boost/type_traits/is_complex.hpp>
#include "bits/PolynomialTermParser.hpp"

namespace PolyEvaluate {
    namespace bmt=boost::math::tools; // for polynomial

    // default template
    template<class T_Roots,class T_Coeff, class Enable = void, class Enable2 = void>
    class PolyEvaluate {
    public:
        explicit PolyEvaluate(const bmt::polynomial<T_Coeff> &poly) {

        }
        T_Roots evalPoly(T_Roots x){
            return (T_Roots)0;
        }

    };

    // specialization for floating point coeffs and roots
    template<class T_Coeff,class T_Roots>
    class PolyEvaluate<T_Coeff,T_Roots,
            typename std::enable_if<(std::is_floating_point<T_Coeff>::value)>::type,
            typename std::enable_if<(std::is_floating_point<T_Roots>::value)>::type> {
    public:
        const bmt::polynomial<T_Coeff>& mPoly;
        explicit PolyEvaluate(const bmt::polynomial<T_Coeff> &poly) : mPoly(poly) {
        }

        /**
         *
         * @param roots
         * @param start
         */
        T_Roots evalPoly(T_Roots x) {
            T_Roots result = mPoly[mPoly.degree()] * x;
            for (int i = (int)(mPoly.degree() - 1); i > 0; i--) {
                result = (result + mPoly[i]) * x;
            }
            result += mPoly[0];
            return result;
        }
    };

    // specialization for floating point coeffs and complex roots
    template<class T_Coeff,class T_Roots>
    class PolyEvaluate<T_Coeff,T_Roots,
            typename std::enable_if<(std::is_floating_point<T_Coeff>::value)>::type,
            typename std::enable_if<(boost::is_complex<T_Roots>::value)>::type> {
    public:
        const bmt::polynomial<T_Coeff>& mPoly;
        explicit PolyEvaluate(const bmt::polynomial<T_Coeff> &poly) : mPoly(poly) {
        }

        /**
         *
         * @param roots
         * @param start
         */
        T_Roots evalPoly(T_Roots x) {
            T_Roots result = (T_Roots)mPoly[mPoly.degree()] * x;
            for (int i = (int)(mPoly.degree() - 1); i > 0; i--) {
                result = (result + (T_Roots)mPoly[i]) * x;
            }
            result += (T_Roots)mPoly[0];
            return  result;
        }
    };
    // specialization for complex coeffs and roots
    template<class T_Coeff,class T_Roots>
    class PolyEvaluate<T_Coeff,T_Roots,
            typename std::enable_if<(boost::is_complex<T_Coeff>::value)>::type,
            typename std::enable_if<(boost::is_complex<T_Roots>::value)>::type> {
    public:
        const bmt::polynomial<T_Coeff>& mPoly;
        explicit PolyEvaluate(const bmt::polynomial<T_Coeff> &poly) : mPoly(poly) {
        }

        /**
         *
         * @param roots
         * @param start
         */
        T_Roots evalPoly(T_Roots x) {
            T_Roots result = (T_Roots)mPoly[mPoly.degree()] * x;
            for (int i = (int)(mPoly.degree() - 1); i > 0; i--) {
                result = (result + (T_Roots)mPoly[i]) * x;
            }
            result += (T_Roots)mPoly[0];
            return  result;
        }
    };
}
#endif //PROYECTO1_POLYEVALUATE_HPP
