//
// Created by chris on 9/3/18.
//

#ifndef PROYECTO1_HELPERMULLER_HPP
#define PROYECTO1_HELPERMULLER_HPP

#include <vector>
#include <complex>
#include <type_traits>
#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>
#include "PolynomialFormulaFormat.hpp"
#include "PolynomialParser.hpp"
#include "PolyEvaluate.hpp"
namespace bmt=boost::math::tools; // for polynomial

namespace helperMuller{
    // default template
    template<class T_Coeff,class T_Roots, class Enable = void, class Enable2 = void>
    class HelperMuller {
    public:
        explicit HelperMuller(const bmt::polynomial<T_Coeff> &poly) {

        }
        T_Roots solve(const T_Roots start = T_Roots()){
            std::cout<< "Nothing" << std::endl;
            return (T_Roots)0;
        }

    };

    // specialization for floating point coeffs and roots
    template<class T_Coeff, class T_Roots>
    class HelperMuller<T_Coeff,T_Roots,
            typename std::enable_if<(std::is_floating_point<T_Coeff>::value)>::type,
            typename std::enable_if<(std::is_floating_point<T_Roots>::value)>::type> {
    public:
        const bmt::polynomial<T_Coeff>& mPoly;
        explicit HelperMuller(const bmt::polynomial<T_Coeff> &poly) : mPoly(poly){

        }

        /**
         *
         * @param roots
         * @param start
         */
        T_Roots solve(const T_Roots start = T_Roots()){

            T_Roots xr,x0,x1,x2,den,diff,dxr;
            int maxit = 1000;

            diff = (0.001);
            xr = start;
            x0 = xr-diff;
            x1 = start;
            x2 = xr+diff;

            for(int i = 2; i < maxit; i++){
                T_Roots h0 = x1 - x0;
                T_Roots h1 = x2 - x1;
                T_Roots d0 = (PolyEvaluate::PolyEvaluate<T_Roots,T_Coeff>(mPoly).evalPoly(x1) - PolyEvaluate::PolyEvaluate<T_Roots,T_Coeff>(mPoly).evalPoly(x0))/h0;
                T_Roots d1 = (PolyEvaluate::PolyEvaluate<T_Roots,T_Coeff>(mPoly).evalPoly(x2) - PolyEvaluate::PolyEvaluate<T_Roots,T_Coeff>(mPoly).evalPoly(x1))/h1;
                T_Roots a = (d1 - d0)/(h1+h0);
                T_Roots b = a*h1 +d1;
                T_Roots c = PolyEvaluate::PolyEvaluate<T_Roots,T_Coeff>(mPoly).evalPoly(x2);
                T_Roots rad = std::sqrt(b*b - ((T_Roots)4.0*a*c));
                T_Roots temp_Pos = std::abs(b+rad);
                T_Roots temp_Neg = std::abs(b-rad);
                if(abs(temp_Pos) > abs(temp_Neg) ){
                    den = b+rad;
                }
                else{
                    den = b-rad;
                }
                dxr = -(T_Roots)2*c/den;

                xr = x2 + dxr;
                if(abs(dxr)<(DBL_EPSILON*abs(xr))){
                    break;
                }
                else{
                    std::cout<< xr  << std::endl;
                    x0 = x1;
                    x1 = x2;
                    x2 = xr;
                }
            }
            std::cout<< "Double roots" << std::endl;
            return xr;
        }
    };

    // specialization for real coefficients with complex roots
    template<class T_Coeff, class T_Roots>
    class HelperMuller<T_Coeff,T_Roots,
            typename std::enable_if<(std::is_floating_point<T_Coeff>::value)>::type,
            typename std::enable_if<(boost::is_complex<T_Roots>::value)>::type> {
    public:
        const bmt::polynomial<T_Coeff>& mPoly;
        explicit HelperMuller(const bmt::polynomial<T_Coeff> &poly) : mPoly(poly){

        }
        T_Roots solve(const T_Roots start = T_Roots()){

            T_Roots xr,x0,x1,x2,den,diff,dxr;
            int maxit = 1000;

            diff = (T_Roots)(0.001);
            xr = start;
            x0 = xr-diff;
            x1 = start;
            x2 = xr+diff;

            for(int i = 2; i < maxit; i++){
                T_Roots h0 = x1 - x0;
                T_Roots h1 = x2 - x1;
                T_Roots d0 = (PolyEvaluate::PolyEvaluate<T_Coeff,T_Roots>(mPoly).evalPoly(x1) - PolyEvaluate::PolyEvaluate<T_Coeff,T_Roots>(mPoly).evalPoly(x0))/h0;
                T_Roots d1 = (PolyEvaluate::PolyEvaluate<T_Coeff,T_Roots>(mPoly).evalPoly(x2) - PolyEvaluate::PolyEvaluate<T_Coeff,T_Roots>(mPoly).evalPoly(x1))/h1;
                T_Roots a = (d1 - d0)/(h1+h0);
                T_Roots b = a*h1 +d1;
                T_Roots c = PolyEvaluate::PolyEvaluate<T_Coeff,T_Roots>(mPoly).evalPoly(x2);
                T_Roots rad = std::sqrt(b*b - (T_Roots)4.0*a*c);
                T_Roots temp_Pos = std::abs(b+rad);
                T_Roots temp_Neg = std::abs(b-rad);
                if(abs(abs(temp_Pos)) > abs(abs(temp_Neg)) ){
                    den = b+rad;
                }
                else{
                    den = b-rad;
                }
                dxr = -(T_Roots)2*c/den;

                xr = x2 + dxr;
                if(abs(abs(dxr))<(DBL_EPSILON*abs(abs(xr)))){
                    break;
                }
                else{
                    std::cout<< xr  << std::endl;
                    x0 = x1;
                    x1 = x2;
                    x2 = xr;
                }
            }
            std::cout<< "Complex Roots" << std::endl;
            return xr;

        }
    };

    // specialization for complex roots and coefficients
    template<class T_Coeff, class T_Roots>
    class HelperMuller<T_Coeff, T_Roots,
            typename std::enable_if<(boost::is_complex<T_Coeff>::value)>::type,
            typename std::enable_if<(boost::is_complex<T_Roots>::value)>::type> {
    public:
        const bmt::polynomial<T_Coeff> &mPoly;

        explicit HelperMuller(const bmt::polynomial<T_Coeff> &poly) : mPoly(poly) {

        }

        T_Roots solve(const T_Roots start = T_Roots()) {
            T_Roots xr,x0,x1,x2,den,diff,dxr;
            int maxit = 1000;

            diff = (T_Roots)(0.001);
            xr = start;
            x0 = xr-diff;
            x1 = start;
            x2 = xr+diff;

            for(int i = 2; i < maxit; i++){
                T_Roots h0 = x1 - x0;
                T_Roots h1 = x2 - x1;
                T_Roots d0 = (PolyEvaluate::PolyEvaluate<T_Coeff,T_Roots>(mPoly).evalPoly(x1) - PolyEvaluate::PolyEvaluate<T_Coeff,T_Roots>(mPoly).evalPoly(x0))/h0;
                T_Roots d1 = (PolyEvaluate::PolyEvaluate<T_Coeff,T_Roots>(mPoly).evalPoly(x2) - PolyEvaluate::PolyEvaluate<T_Coeff,T_Roots>(mPoly).evalPoly(x1))/h1;
                T_Roots a = (d1 - d0)/(h1+h0);
                T_Roots b = a*h1 +d1;
                T_Roots c = PolyEvaluate::PolyEvaluate<T_Coeff,T_Roots>(mPoly).evalPoly(x2);
                T_Roots rad = std::sqrt(b*b - (T_Roots)4.0*a*c);
                T_Roots temp_Pos = std::abs(b+rad);
                T_Roots temp_Neg = std::abs(b-rad);
                if(abs(temp_Pos) > abs(temp_Neg)){
                    den = b+rad;
                }
                else{
                    den = b-rad;
                }
                dxr = -(T_Roots)2*c/den;

                xr = x2 + dxr;
                if(abs(dxr)<(DBL_EPSILON*abs(xr))){
                    break;
                }
                else{
                    std::cout<< xr  << std::endl;
                    x0 = x1;
                    x1 = x2;
                    x2 = xr;
                }
            }
            std::cout<< "Complex roots end coeff" << std::endl;
            return xr;
        }
    };
}


#endif //PROYECTO1_HELPERMULLER_HPP