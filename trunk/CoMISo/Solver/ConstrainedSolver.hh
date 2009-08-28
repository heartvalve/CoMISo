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
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  CoMISo is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with CoMISo.  If not,                                      *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/ 


//=============================================================================
//
//  CLASS ConstrainedSolver
//
//=============================================================================


#ifndef ACG_CONSTRAINEDSOLVER_HH
#define ACG_CONSTRAINEDSOLVER_HH


//== INCLUDES =================================================================
#include <CoMISo/Config/CoMISoDefines.hh>

#include "GMM_Tools.hh"
#include "MISolver.hh"
#include <vector>

//== FORWARDDECLARATIONS ======================================================

//== DEFINES ==================================================================
#define ROUND(x) ((x)<0?int((x)-0.5):int((x)+0.5))
//== NAMESPACES ===============================================================

namespace ACG {
//== CLASS DEFINITION =========================================================

/** \class ConstrainedSolver ConstrainedSolver.hh <ACG/.../ConstrainedSolver.hh>

  Takes a linear (symmetric) system of equations and a set of linear constraints and solves it.
 */

class COMISODLLEXPORT ConstrainedSolver
{
public:
  typedef gmm::csc_matrix<double>  CSCMatrix;


  /// default Constructor
  ConstrainedSolver()  { }

  /// Destructor
  ~ConstrainedSolver() { }

/** @name Contrained solvers
 * Functions to solve constrained linear systems of the form Ax=b (stemming from quadratic energies). 
 * The constraints can be linear constraints of the form \f$ x_1*c_1+ \cdots +x_n*c_n=c \f$ as well as integer constraints \f$x_i\in \mathbf{Z}\f$. 
 * The system is then solved with these constraints in mind. For solving the system the Mixed-Integer Solver \a MISolver is used. 
 */
/*@{*/

/// Quadratic matrix constrained solver
/**  
  *  Takes a system of the form Ax=b, a constraint matrix C and a set of variables _to_round to be rounded to integers. \f$ A\in \mathbf{R}^{n\times n}\f$
  *  @param _constraints row matrix with rows of the form \f$ [c_1, c_2, \cdots, c_n, c_{n+1}] \f$ corresponding to the linear equation \f$ c_1*x_1+\cdots+c_n*x_n + c_{n+1}=0 \f$.
  *  @param _A nxn-dimensional column matrix of the system 
  *  @param _x n-dimensional variable vector
  *  @param _rhs n-dimensional right hand side.
  *  @param _idx_to_round indices i of variables x_i that shall be rounded
  *  @param _reg_factor regularization factor. Helps unstable, low rank system to become solvable. Adds \f$ \_reg\_factor*mean(trace(_A))*Id \f$ to A.
  *  @param _show_miso_settings should the (QT) dialog of the Mixed-Integer solver pop up?
  *  @param _show_timings shall some timings be printed?
  */
  template<class RMatrixT, class CMatrixT, class VectorT, class VectorIT >
  void solve(
      RMatrixT& _constraints,
      CMatrixT& _A, 
      VectorT&  _x,
      VectorT&  _rhs,
      VectorIT& _idx_to_round,
      double    _reg_factor = 0.0,
      bool      _show_miso_settings = true,
      bool      _show_timings = true );


/// Non-Quadratic matrix constrained solver
/**  
  *  Same as above, but performs the elimination of the constraints directly on the B matrix of \f$ x^\top B^\top Bx \f$, where B has m rows (equations) and (n+1) columns \f$ [x_1, x_2, \dots, \x_n, -rhs] \f$. \note This function might be more efficient in some cases, but generally the solver for the quadratic matrix above is a safer bet. Needs further testing.
  *  \note Internally the \f$ A=B^\top B \f$ matrix is formed.
  *  @param _constraints row matrix with rows of the form \f$ [c_1, c_2, \cdots, c_n, c_{n+1}] \f$ corresponding to the linear equation \f$ c_1*x_1+\cdots+c_n*x_n + c_{n+1}=0 \f$.
  *  @param _B mx(n+1)-dimensional column matrix of the system 
  *  @param _x n-dimensional variable vector
  *  @param _idx_to_round indices i of variables x_i that shall be rounded
  *  @param _reg_factor regularization factor. Helps unstable, low rank system to become solvable.
  *  @param _show_miso_settings should the (QT) dialog of the Mixed-Integer solver pop up?
  *  @param _show_timings shall some timings be printed?
  */
  template<class RMatrixT, class VectorT, class VectorIT >
  void solve(
      RMatrixT& _constraints,
      RMatrixT& _B, 
      VectorT&  _x,
      VectorIT& _idx_to_round,
      double    _reg_factor = 0.0,
      bool      _show_miso_settings = true,
      bool      _show_timings = true );
/*@}*/

/** @name Eliminate constraints
 * Functions to eliminate (or integrate) linear constraints from an equation system. These functions are used internally by the \a solve functions.
 */
/*@{*/

/// Make constraints independent
/**  
  *  This function performs a Gauss elimination on the constraint matrix making the constraints easier to eliminate. 
  *  \note A certain amount of independence of the constraints is assumed.
  *  \note contradicting constraints will be ignored.
  *  \warning care must be taken when non-trivial constraints occur where some of the variables contain integer-variables (to be rounded) as the optimal result might not always occur.
  *  @param _constraints  row matrix with constraints
  *  @param _idx_to_round indices of variables to be rounded (these must be considered.)
  *  @param _c_elim the "returned" vector of variable indices and the order in which the can be eliminated.
  */
  template<class RMatrixT, class VectorIT >
  void make_constraints_independent(
      RMatrixT&         _constraints,
			VectorIT&         _idx_to_round,
			std::vector<int>& _c_elim );

/// Eliminate constraints on a factored matrix B
/**  
  *  \note Constraints are assumed to have been made independent by \a make_constraints_independent.
  *  @param _constraints row matrix with constraints (n+1 columns) 
  *  @param _B system row matrix mx(n+1)
  *  @param _idx_to_round indices to be rounded
  *  @param _c_elim the indices of the variables to be eliminated.
  *  @param _new_idx the created re-indexing map. new_idx[i] = -1 means x_i eliminated, new_idx[i] = j means x_i is now at index j.
  *  @param _Bcol resulting (smaller) column matrix to be used for future computations. (e.g. convert to CSC and solve)
  */
  template<class SVector1T, class SVector2T, class VectorIT, class SVector3T>
  void eliminate_constraints(
      gmm::row_matrix<SVector1T>& _constraints,
			gmm::row_matrix<SVector2T>& _B, 
			VectorIT&                   _idx_to_round,
			std::vector<int>&           _c_elim,
			std::vector<int>&           _new_idx,
			gmm::col_matrix<SVector3T>& _Bcol);

/// Eliminate constraints on a quadratic matrix A
/**  
  *  \note Constraints are assumed to have been made independent by \a make_constraints_independent.
  *  \note _x must have correct size (same as _rhs)
  *  @param _constraints row matrix with constraints (n+1 columns) 
  *  @param _A system row matrix nxn)
  *  @param _x variable vector
  *  @param _rhs right hand side
  *  @param _idx_to_round indices to be rounded
  *  @param _c_elim the indices of the variables to be eliminated.
  *  @param _new_idx the created re-indexing map. new_idx[i] = -1 means x_i eliminated, new_idx[i] = j means x_i is now at index j.
  *  @param _Acsc resulting (smaller) column (csc) matrix to be used for future computations.
  */
 
  template<class SVector1T, class SVector2T, class VectorIT, class CSCMatrixT>
  void eliminate_constraints(
      gmm::row_matrix<SVector1T>& _constraints,
      gmm::col_matrix<SVector2T>& _A, 
      std::vector<double>&        _x, 
      std::vector<double>&        _rhs, 
      VectorIT&                   _idx_to_round,
      std::vector<int>&           _c_elim,
      std::vector<int>&           _new_idx,
      CSCMatrixT&                 _Acsc);

/// Restore a solution vector to the un-eliminated size
/**  
  *  @param _constraints row matrix with constraints (n+1 columns) 
  *  @param _x solution vector to reduced/eliminated system (result will also be written here)
  *  @param _c_elim vector of eliminated indices 
  *  @param _new_idx re-indexing vector
  */
 
  template<class RMatrixT, class VectorT >
  void restore_eliminated_vars( RMatrixT&         _constraints,
				VectorT&          _x,
				std::vector<int>& _c_elim,
				std::vector<int>& _new_idx);


/*@}*/


/** @name Verify the result.
 * Functions to verify the result of the constrained solver. Are the constraints met, are the correct variables correctly rounded ...
 */
/*@{*/


  template<class RMatrixT, class CMatrixT, class VectorT>
  double verify_constrained_system( 
				   const RMatrixT& _conditions,
				   const CMatrixT& _A,
				   const VectorT&  _x,
				   const VectorT&  _rhs);

  template<class RMatrixT, class CMatrixT, class VectorT, class VectorIT>
  double verify_constrained_system_round( 
					 const RMatrixT& _conditions,
					 const CMatrixT& _A,
					 const VectorT&  _x,
					 const VectorT&  _rhs,
					 const VectorIT& _idx_to_round);

  template<class RMatrixT, class VectorT, class VectorIT>
  void verify_mi_factored( const RMatrixT& _conditions,
			   const RMatrixT& _B, 
			   const VectorT&  _x,
			   const VectorIT& _idx_to_round );
/*@}*/

private:

  template<class RowT, class MatrixT>
  void add_row( int       _row_i,
		double    _coeff,
		RowT      _row, 
		MatrixT&  _mat );

  template<class RowT, class RMatrixT, class CMatrixT>
  void add_row_simultaneously( int       _row_i,
			       double    _coeff,
			       RowT      _row, 
			       RMatrixT& _rmat,
			       CMatrixT& _cmat );


  template<class CMatrixT, class VectorT, class VectorIT>
  double setup_and_solve_system( CMatrixT& _B,
			       VectorT&  _x,
			       VectorIT& _idx_to_round,
			       double    _reg_factor,
			       bool      _show_miso_settings);


  // warning: order of replacement not the same as in _columns (internal sort)
  template<class CMatrixT>
  void eliminate_columns( CMatrixT& _M,
			  const std::vector< int >& _columns);

private:

  /// Copy constructor (not used)
  ConstrainedSolver(const ConstrainedSolver& _rhs);

  /// Assignment operator (not used)
  ConstrainedSolver& operator=(const ConstrainedSolver& _rhs);
};


//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_CONSTRAINEDSOLVER_C)
#define ACG_CONSTRAINEDSOLVER_TEMPLATES
#include "ConstrainedSolverT.cc"
#endif
//=============================================================================
#endif // ACG_CONSTRAINEDSOLVER_HH defined
//=============================================================================

