//
// Created by allan on 02/09/18.
//

#ifndef HELPER_INSERT_H
#define HELPER_INSERT_H
#include "Rpoly.hpp"
#include <vector>
#include <complex>
#include <type_traits>
#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>

namespace bmt=boost::math::tools; // for polynomial

namespace hvins{
    // default template
    template<class T_Coeff,class T_Roots, class Enable = void, class Enable2 = void>
    class Helper {
    public:
        explicit Helper(const bmt::polynomial<T_Coeff> &poly) {

        }
        void test(std::vector<T_Roots>  &roots){

        }
    };

    // specialization for floating point coeffs and roots
    template<class T_Coeff, class T_Roots>
    class Helper<T_Coeff,T_Roots,
            typename std::enable_if<(std::is_floating_point<T_Coeff>::value)>::type,
            typename std::enable_if<(std::is_floating_point<T_Roots>::value)>::type> {
    public:
        const bmt::polynomial<T_Coeff>& mPoly;
        explicit Helper(const bmt::polynomial<T_Coeff> &poly) : mPoly(poly){

        }
        void test(std::vector<T_Roots>  &roots){
            unsigned long degree = mPoly.degree();
            std::vector<T_Coeff> zeror(degree), zeroi(degree);
            Jenkins::RPoly<T_Coeff>().findRoots(mPoly.data(), static_cast<int>(degree), zeror, zeroi);
            for (unsigned long i = 0; i < degree; ++i) {
                roots.push_back(zeror[i]);
            }
        }
    };

    // specialization for complex roots
    template<class T_Coeff, class T_Roots>
    class Helper<T_Coeff,T_Roots,
            typename std::enable_if<(std::is_floating_point<T_Coeff>::value)>::type,
            typename std::enable_if<(boost::is_complex<T_Roots>::value)>::type> {
    public:
        const bmt::polynomial<T_Coeff>& mPoly;
        explicit Helper(const bmt::polynomial<T_Coeff> &poly) : mPoly(poly){

        }
        void test(std::vector<T_Roots>  &roots){
            unsigned long degree = mPoly.degree();
            std::vector<T_Coeff> zeror(degree), zeroi(degree);
            Jenkins::RPoly<T_Coeff>().findRoots(mPoly.data(), static_cast<int>(degree), zeror, zeroi);
            for (unsigned long i = 0; i < degree; ++i) {
                roots.push_back(T_Roots(zeror[i],zeroi[i]));
            }
        }
    };
}


#endif //HELPER_INSERT_H
