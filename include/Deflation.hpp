/**
 * Copyright (C) 2017-2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @author Pablo Alvarado
 * @date   18.08.2018
 */

#ifndef ANPI_DEFLATION_HPP
#define ANPI_DEFLATION_HPP

#include <vector>
#include <type_traits>

#include <boost/type_traits/is_complex.hpp>
#include <boost/math/tools/polynomial.hpp>


namespace anpi {
    namespace bmt=boost::math::tools; // bt as alias for boost::math::tools

    /**
     * Deflate polynomial
     *
     * @param[in] poly Input polynomial to be deflated with the provided root
     * @param[in] root Root of poly to be deflated with.
     * @param[out] residuo Residual of the polynomial deflation
     * @return deflated polynomial
     */
    template<class T>
    bmt::polynomial<T> deflate(const bmt::polynomial<T> &poly,
                               const T &root,
                               T &residue
                               ) {

        /*
         * Important data is defined, the degree of the function is
         * stored as a constant, the coefficients of the polynomial fuction
         * are saved and a new vector to store the results is generated
         */
        int degree = poly.degree();
        std::vector<T> term = poly.data();
        std::vector<T> factorTerms(degree);


        //Inicialization of the variables to use in the deflation process
        factorTerms[degree] = 0;
        residue = term[degree];

        //Deflation cycle
        for (int i = degree - 1; i >= 0; --i) {

            factorTerms[i] = residue;             // The coefficient is stored
            residue = term[i] + residue * root;     // The new coefficient is calculated

        }

        // Store the output in polynomial
        bmt::polynomial<T> const factoredPoly(factorTerms.begin(), factorTerms.end());

        return factoredPoly;


    }

    /**
     * Deflate polynomial with a second order polynomial.
     *
     * The second order polynomial equals x^2 -2 Re(root)x + |root|^2.
     *
     * @param[in] poly Input polynomial to be deflated with the provided root
     * @param[in] root Root of poly to be deflated with.
     * @param[out] residuo Residual of the polynomial deflation
     * @return deflated polynomial
     */
    template<class T>
    bmt::polynomial<T> deflate2(const bmt::polynomial<T> &poly,
                                const T &root,
                                T &residue) {

        /*
         * Important data is defined, the degree of the function is
         * stored as a constant, the coefficients of the polynomial fuction
         * are saved and a new vector to store the results is generated.
         *
         * The conjugated root is also saved
         */
        int degree = poly.degree();
        std::complex<T> rootConj = std::conj(root);
        std::vector<T> term = poly.data();
        std::vector<T> factorTerms(degree);


        // Inicialization of the variables to use in the deflation process
        residue = term[degree];
        T factor = T(0);

        // Deflation cycle
        for (int i = degree - 1; i >= 0; --i) {

            factorTerms[i] = factor;                      // The coefficient is stored
            factor = (residue + factor * rootConj).real();  // The new coefficient is stored
            residue = term[i] + residue * root;

        }


        bmt::polynomial<T> const factoredPoly(factorTerms.begin(), factorTerms.end());
        return factoredPoly;


    }


}



#endif
