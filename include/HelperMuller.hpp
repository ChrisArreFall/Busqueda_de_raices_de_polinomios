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
#include "Deflation.hpp"
namespace bmt=boost::math::tools; // for polynomial

namespace helperMuller{
    // default template
    template<class T_Coeff,class T_Roots, class Enable = void, class Enable2 = void>
    class HelperMuller {
    public:
        explicit HelperMuller(const bmt::polynomial<T_Coeff> &poly) {

        }
        T_Roots solve(const T_Roots start = T_Roots(),bmt::polynomial<T_Coeff> &poly = T_Coeff()){
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
        T_Roots solve(const T_Roots start = T_Roots(),bmt::polynomial<T_Coeff> &poly = T_Coeff()){
            //Muller algorithm
            std::complex<T_Roots> xr,x0,x1,x2,den,diff,dxr;
            int maxit = 1000;

            diff = (0.001);
            xr = (T_Roots)start;
            x0 = xr-diff;
            x1 = (T_Roots)start;
            x2 = xr+diff;

            for(int i = 2; i < maxit; i++){
                std::complex<T_Roots> h0 = x1 - x0;
                std::complex<T_Roots> h1 = x2 - x1;
                std::complex<T_Roots> d0 = (PolyEvaluate::PolyEvaluate<std::complex<T_Coeff>,std::complex<T_Roots>>(mPoly).evalPoly(x1) - PolyEvaluate::PolyEvaluate<std::complex<T_Coeff>,std::complex<T_Roots>>(mPoly).evalPoly(x0))/h0;
                std::complex<T_Roots> d1 = (PolyEvaluate::PolyEvaluate<std::complex<T_Coeff>,std::complex<T_Roots>>(mPoly).evalPoly(x2) - PolyEvaluate::PolyEvaluate<std::complex<T_Coeff>,std::complex<T_Roots>>(mPoly).evalPoly(x1))/h1;
                std::complex<T_Roots> a = (d1 - d0)/(h1+h0);
                std::complex<T_Roots> b = a*h1 +d1;
                std::complex<T_Roots> c = PolyEvaluate::PolyEvaluate<std::complex<T_Coeff>,std::complex<T_Roots>>(mPoly).evalPoly(x2);
                std::complex<T_Roots> rad = std::sqrt(b*b - ((T_Roots)4.0*a*c));
                std::complex<T_Roots> temp_Pos = std::abs(b+rad);
                std::complex<T_Roots> temp_Neg = std::abs(b-rad);
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
                    x0 = x1;
                    x1 = x2;
                    x2 = xr;
                }
            }
            //------------------------------------------------------------
            //Deflacion
            T_Coeff tempRes = T_Coeff(0);
            bmt::polynomial<T_Coeff> factoredPoly = anpi::deflate(mPoly,(T_Coeff)xr.real(),tempRes);
            //----------------------------------------------------------------------------------
            poly = factoredPoly;

            return xr.real();
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
        T_Roots solve(const T_Roots start = T_Roots(),bmt::polynomial<T_Coeff> &poly = T_Coeff()){
            //Muller algorithm
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
                    x0 = x1;
                    x1 = x2;
                    x2 = xr;
                }
            }
            //------------------------------------------------------------
            //Deflacion
            T_Coeff residue;

            int degree = mPoly.degree();
            std::vector<T_Coeff> term = mPoly.data();
            std::vector<T_Coeff> factorTerms(degree);


            //Inicialization of the variables to use in the deflation process
            factorTerms[degree] = 0;
            residue = term[degree];

            //Deflation cycle
            for (int i = degree - 1; i >= 0; --i) {

                factorTerms[i] = residue;             // The coefficient is stored
                residue = term[i] + residue * xr.real();     // The new coefficient is calculated

            }

            // Store the output in polynomial
            bmt::polynomial<T_Coeff> const factoredPoly(factorTerms.begin(), factorTerms.end());

            //----------------------------------------------------------------------------------
            poly = factoredPoly;
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

        T_Roots solve(const T_Roots start = T_Roots(),bmt::polynomial<T_Coeff> &poly = T_Coeff()) {
            //Muller algorithm
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
                    x0 = x1;
                    x1 = x2;
                    x2 = xr;
                }
            }
            //------------------------------------------------------------
            //Deflacion
            T_Coeff xrTemp = T_Coeff(xr.real(),xr.imag());
            T_Coeff tempRes = T_Coeff(0,0);
            bmt::polynomial<T_Coeff> factoredPoly = anpi::deflate(mPoly,xrTemp,tempRes);
            //----------------------------------------------------------------------------------
            poly = factoredPoly;
            return xr;
        }
    };
}


#endif //PROYECTO1_HELPERMULLER_HPP
