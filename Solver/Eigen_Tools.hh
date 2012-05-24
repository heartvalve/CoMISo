/*===========================================================================*\
 *                                                                           *
 *                               CoMISo                                      *
 *      Copyright (C) 2008-2009 by Computer Graphics Group, RWTH Aachen      *
 *                           www.rwth-graphics.de                            *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *  This file is part of CoMISo.                                             *
 *                                                                           *
 *  CoMISo is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  CoMISo is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with CoMISo.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                           *
\*===========================================================================*/ 


#ifndef COMISO_Eigen_TOOLS_HH
#define COMISO_Eigen_TOOLS_HH


//== INCLUDES =================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

//#ifdef COMISO_Eigen3_AVAILABLE
//#include <Eigen/Eigen>
//#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
//#include <Eigen/Sparse>
//#endif

#ifndef COMISO_NCHOLMOD
#include <cholmod.h>
#endif


//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO_EIGEN
{

/** \class EigenTools Eigen_Tools.hh

    A collection of helper functions for manipulating (Eigen) matrices.
*/


//== FUNCTION DEFINITION ======================================================

/// Get matrix data (CSC matrix format) from matrix
/** Used by Cholmod wrapper  
 *  @param _mat matrix
 *  @param _c uplo parameter (l, L, u, U, c, C)
 *  @param _values values vector
 *  @param _rowind row indices 
 *  @param _colptr column pointer  */
template<class MatrixT, class REALT, class INTT>
void get_ccs_symmetric_data( const MatrixT&      _mat,
                             const char          _c,
                             std::vector<REALT>& _values,
                             std::vector<INTT>&  _rowind,
                             std::vector<INTT>&  _colptr );

/// Inspect the matrix (print)
/** Prints useful matrix informations such as, dimension, symmetry, zero_rows, zero_cols, nnz, max, min, max_abs, min_abs, NAN, INF
  * @param _A matrix */
template<class MatrixT>
void inspect_matrix( const MatrixT& _A);

/** checks for symmetry
  * @param _A matrix 
  * @return symmetric? (bool)*/
template<class MatrixT>
bool is_symmetric( const MatrixT& _A);



//=============================================================================
} // namespace COMISO_Eigen
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(COMISO_Eigen_TOOLS_C)
#define COMISO_Eigen_TOOLS_TEMPLATES
#include "Eigen_Tools.cc"
#endif
//=============================================================================
#endif // Eigen_TOOLS_HH defined
//=============================================================================

