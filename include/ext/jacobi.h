
// Copyright (C) 2006 Paulo José Saiz Jabardo
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

/* Author:  Paulo J. Saiz Jabardo */

/** \file jacobi.h
  \brief Header file for jacobi library
  
  The jacobi library implements Jacobi polynomials and Gauss-Jacobi quadrature related functions.
  This library was developed to be used in high order finite element methods (Spectral Element/hp
  methods. The interfaces were inspired by the software Polylib (http://www.nektar.info).

  The main sources used to develop this library where:

  (1) Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
      Dover Publications 1965.

  (2) Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
      Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.

  This software is supposed to be an add-on module for the GNU Scientific Library
*/

#ifndef __JACOBI_H__
#define __JACOBI_H__

  /** 
   * \brief Enumeration gor differing types of Gauss quadrature.
   * The gauss_quad_type is used to determine the boundary condition
   * modifications applied to orthogonal polynomials for quadrature
   * rules.
   */
  enum gauss_quad_type
  {
    Gauss,             ///< Gauss quadrature
    Gauss_Lobatto,     ///< Gauss-Lobatto quadrature
    Gauss_Radau_lower, ///< Gauss-Radau quadrature including the node -1
    Gauss_Radau_upper  ///< Gauss-Radau quadrature including the node +1
  };


/**
 * \brief Struct used to store information and memory on the quadrature
 *
 * This strucuture is used to store data about the quadrature. It can be later used to integrate, 
 * derive or interpolate functions. It should be created using the jac_quadrature_alloc function.
 */
template<typename Tp>
  struct jac_quadrature
  {
    /// Quadrature type
    enum gauss_quad_type type;

    /// Number of quadrature points
    int Q;

    /// Alpha weight of the quadrature
    Tp alpha;

    /// Beta weight of the quadrature
    Tp beta;

    jac_quadrature(gauss_quad_type qtype, int nq, Tp a, Tp b)
    : type(qtype), Q(nq), alpha(a), beta(b),
      x(nq), w(nq), D(nq),
      xp{}, imat{}
    {
      this->quadrature_zwd();
    }

    /// Calculates the integral of a function with values at quadrature points
    /// given by f with quadrature weights w.
    template<typename Func>
      Tp integrate(Func fun);

    ~jac_quadrature() = default;

  private:

    /// Calculates the quadrature zeros, weights and derivative matrix.
    int quadrature_zwd();

    /// Allocates memory for an interpolation matrix.
    int interpmat_alloc(int npoints, Tp* xp);

    /// Calculates the derivative at quadrature points of a function f
    /// and derivative matrix D
    int differentiate(Tp* f, Tp* d);

    /// Interpolates the function given by f using the interpolation matrix.
    int interpolate(Tp* f, Tp* fout);

    /// Array that stores the nodes coordinates
    std::vector<Tp> x;

    /// Array that stores quadrature weights
    std::vector<Tp> w;

    /// Array that stores the derivative matrix
    std::vector<Tp> D;

    /// Interpolation points
    std::vector<Tp> xp;

    /// Array That stores the interpolation matrix
    std::vector<Tp> imat;


    /// Calculates the derivative of the Jacobi polynomial of order n
    Tp jacobi_deriv(Tp x, int n, Tp alpha, Tp beta);

    /// Calculates Jacobi polynomials at an array of points
    //int jacobi_value_array(int nx, Tp* x, Tp* result_array, int n, Tp alpha, Tp beta);

    /// Calculates the derivative of Jacobi polynomials at an array of points
    //int jacobi_deriv_array(int nx, Tp* x, Tp* result_array, int n, Tp alpha, Tp beta);

    /// Calculates the zeros of Jacobi polynomials in the interval -1 up to 1
    int jacobi_zeros(Tp* x, int m, Tp alpha, Tp beta);


    /// Calculates the Gauss-Jacobi quadrature points.
    int zeros_gj();

    /// Calculates the Gauss-Jacobi quadrature weights.
    int weights_gj();

    /// Calculates the Gauss-Jacobi derivative matrix.
    int diffmat_gj();

    /// Calculates the Lagrange polynomials through Gauss-Jacobi quadratures.
    Tp lagrange_gj(int i, Tp zz);

    /// Calculates the interpolation matrix for GJ quadrature.
    int interpmat_gj();

    // **** GAUSS-LOBATTO QUADRATURE ****

    /// Calculates the Gauss-Lobatto-Jacobi quadrature points.
    int zeros_glj();

    /// Calculates the Gauss-Lobatto-Jacobi quadrature weights.
    int weights_glj();

    /// Calculates the Gauss-Lobatto-Jacobi derivative matrix.
    int diffmat_glj();

    /// Calculates the Lagrange polynomials through Gauss-Lobatto-Jacobi quadratures
    Tp lagrange_glj(int i, Tp zz);

    /// Calculates the interpolation matrix for GLJ quadrature
    int interpmat_glj();

    // **** GAUSS-RADAU -1 QUADRATURE ****

    /// Calculates the Gauss-Radau-Jacobi quadrature points (point -1 included).
    int zeros_grjm();

    /// Calculates the Gauss-Radau-Jacobi quadrature weights (point -1 included).
    int weights_grjm();

    /// Calculates the Gauss-Radau-Jacobi derivative matrix (point -1 included).
    int diffmat_grjm();

    /// Calculates the Lagrange polynomials through Gauss-Radau-Jacobi quadratures
    /// (point -1 included)
    Tp lagrange_grjm(int i, Tp zz);

    /// Calculates the interpolation matrix for GRJM quadrature.
    int interpmat_grjm();

    // **** GAUSS-RADAU +1 QUADRATURE ****

    /// Calculates the Gauss-Radau-Jacobi quadrature points (point +1 included).
    int zeros_grjp();

    /// Calculates the Gauss-Radau-Jacobi quadrature weights (point +1 included).
    int weights_grjp();

    /// Calculates the Gauss-Radau-Jacobi derivative matrix (point +1 included).
    int diffmat_grjp();

    /// Calculates the Lagrange polynomials through Gauss-Radau-Jacobi quadratures
    /// (point +1 included).
    Tp lagrange_grjp(int i, Tp zz);

    /// Calculates the interpolation matrix for GRJP quadrature.
    int interpmat_grjp();
  };

#include <ext/gauss_jacobi_interface.tcc>
#include <ext/gauss_jacobi_integrate.tcc>

#endif /* __JACOBI_H__ */
