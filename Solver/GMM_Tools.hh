//=============================================================================
//
//  CLASS GMM_Tools
//  Author:  David Bommes
//           Henrik Zimmer
//  WWW:     www.rwth-graphics.de
//
// Version: $Revision: 1$
// Date:    $Date: 27-05-2009$
//
//=============================================================================


#ifndef GMM_GMM_TOOLS_HH
#define GMM_GMM_TOOLS_HH


//== INCLUDES =================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <gmm/gmm.h>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace gmm
{

/** \class GMMTools GMM_Tools.hh

    A collection of helper functions for manipulating (gmm) matrices.
*/


//== FUNCTION DEFINITION ======================================================

/** @name Variable elimination
 * These functions are used to eliminate (one or more) variables x_i from an equation system Ax=b. Elimination meaning that x_i has been assigned a value x_i = c and is considered a constant, this changes entries of the matrix which depend on i and finally "eliminates" the ith row and column and updates the rhs. */
/*@{*/

/// Eliminate multiple variables from a CSC matrix.
/**  
 *  \note eliminate_csc_vars2 is probably more efficient
 *  @param _evar indices of variables to be eliminated
 *  @param _eval values c_i of x_i to be eliminated, x_i = c_i
 *  @param _A CSC Matrix of the equation system
 *  @param _x variable vector of equation system
 *  @param _rhs right-hand side vector of equation system  */
template<class ScalarT, class VectorT, class RealT, class IntegerT>
void eliminate_csc_vars(
    const std::vector<IntegerT>&     _evar,
    const std::vector<ScalarT>&      _eval,
    typename gmm::csc_matrix<RealT>&  _A,
    VectorT&                          _x,
    VectorT&                          _rhs );

/// Eliminate variables from a CSC matric.
template<class ScalarT, class VectorT, class RealT, class IntegerT>
void eliminate_csc_vars2(
    const std::vector<IntegerT>&     _evar,
    const std::vector<ScalarT>&      _eval,
    typename gmm::csc_matrix<RealT>&  _A,
    VectorT&                          _x,
    VectorT&                          _rhs );

/// Eliminate only one variable x_i = c (CSC matrices)
/** Specialization to eliminating one varaible
 *  @param _j index of variable to be eliminated
 *  @param _value_j value c of x_i to be eliminated, x_i = c
 *  @param _A CSC Matrix of the equation system
 *  @param _x variable vector of equation system
 *  @param _rhs right-hand side vector of equation system */
template<class ScalarT, class VectorT, class RealT>
void eliminate_var( 
    const unsigned int                _j,
    const ScalarT                     _value_j,
    typename gmm::csc_matrix<RealT>&  _A,
    VectorT&                          _x,
    VectorT&                          _rhs );



/// eliminate multiple variables from a (NON CSC) linear system by fixin x[j] = _value_j
/**  
 *  @param _evar indices of variables to be eliminated
 *  @param _eval values c_i of x_i to be eliminated, x_i = c_i
 *  @param _A (non-CSC) Matrix of the equation system
 *  @param _x variable vector of equation system
 *  @param _rhs right-hand side vector of equation system */
template<class IntegerT, class ScalarT, class VectorT, class MatrixT>
void eliminate_vars(
    const std::vector<IntegerT>& _evar,
    const std::vector<ScalarT>&  _eval,
    MatrixT&                     _A,
    VectorT&                     _x,
    VectorT&                     _rhs );

/// Eliminate only one variable x_i = c (non-CSC matrices)
/** Specialization to eliminating one varaible
 *  @param _j index of variable to be eliminated
 *  @param _value_j value c of x_i to be eliminated, x_i = c
 *  @param _A (non-CSC) Matrix of the equation system
 *  @param _x variable vector of equation system
 *  @param _rhs right-hand side vector of equation system */
template<class ScalarT, class VectorT, class SVT>
void eliminate_var(
    const unsigned int     _j,
    const ScalarT          _value_j,
    gmm::col_matrix<SVT>&  _A,
    VectorT&               _x,
    VectorT&               _rhs );

/// update indices of eliminated variables 
/** 
 *  @param _evar indices of variables that were eliminated
 *  @param _idx index map to be changed.
 *  @param _dummy value to which eliminated entries of _idx are set.  */
template<class IntegerT, class IntegerT2>
void eliminate_vars_idx( 
    const std::vector<IntegerT >&  _evar,
    std::vector<IntegerT2>&        _idx,
    IntegerT2                      _dummy = -1 );

/// update index of eliminated variable
/**  Specialization to update only one eliminated variable
 *  @param _evar index of variable that was eliminated
 *  @param _idx index map to be changed.
 *  @param _dummy value to which eliminated entry of _idx is set.  */
template<class IntegerT, class IntegerT2>
void eliminate_var_idx( 
    const IntegerT          _evar,
    std::vector<IntegerT2>& _idx,
    IntegerT2               _dummy = -1 );

/*@}*/


/// Get matrix data (CSC matrix format) from matrix
/** Used by Cholmod wrapper  
 *  @param _mat matrix
 *  @param _c uplo parameter (l, L, u, U)
 *  @param _values values vector
 *  @param _rowind row indices 
 *  @param _colptr column pointer  */
template<class MatrixT, class REALT, class INTT>
void get_ccs_symmetric_data( const MatrixT&      _mat,
                             const char          _c,
                             std::vector<REALT>& _values,
                             std::vector<INTT>&  _rowind,
                             std::vector<INTT>&  _colptr );

/// Regularize matrix
/**  Makes matrices with rank(_mat)<n solvable. 
  *  Add factor*avg(trace(_mat))*Identity to _mat.
 *  @param _mat Matrix to regularize
 *  @param _v factor in factor*avg(trace(_mat))*Identity  */
template<class MatrixT>
void regularize_hack( MatrixT& _mat, double _v = 1e-6 );


/// Local Gauß Seidel update of lin. equation system.
/**  
 *  Add factor*avg(trace(_mat))*Identity to _mat.
 *  @param _A Matrix of linear system
 *  @param _x variable vector of linear system
 *  @param _rhs right hand side of linear system
 *  @param _max_iter maximum number of iterations
 *  @param _tolerance error tolerance threshold
 *  @return number of iterations performed */
template<class MatrixT, class VectorT>
int gauss_seidel_local(
    MatrixT&                  _A,
    VectorT&                  _x,
    VectorT&                  _rhs, 
    std::vector<unsigned int> _idxs,
    int                       _max_iter = 10000,
    double                    _tolerance = 1e-6 );


/// Residuum norm of linear system
/** residuum = Ax-b
  * @param _A Matrix
  * @param _x Variables
  * @param _rhs right hand side
  * @return norm Ax-rhs */
template<class MatrixT, class VectorT>
double residuum_norm( MatrixT& _A, VectorT& _x, VectorT& _rhs );


/// Inspect the matrix (print)
/** Prints useful matrix informations such as, dimension, symmetry, zero_rows, zero_cols, nnz, max, min, max_abs, min_abs, NAN, INF
  * @param _A matrix */
template<class MatrixT>
void inspect_matrix( const MatrixT& _A);



//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(GMM_GMM_TOOLS_C)
#define GMM_GMM_TOOLS_TEMPLATES
#include "GMM_Tools.cc"
#endif
//=============================================================================
#endif // GMM_GMM_TOOLS_HH defined
//=============================================================================

