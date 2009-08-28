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
//  CLASS MISolver
//
//=============================================================================


#ifndef ACG_MISOLVER_HH
#define ACG_MISOLVER_HH


//== INCLUDES =================================================================
#include <CoMISo/Config/CoMISoDefines.hh>

#include "GMM_Tools.hh"
#include "CholmodSolver.hh"

#include <vector>

#define ROUND_MI(x) ((x)<0?int((x)-0.5):int((x)+0.5))


//== FORWARDDECLARATIONS ======================================================


namespace ACG {
class MISolverDialog;
}

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================



/** \class MISolver MISolver.hh 

    Mixed-Integer Solver.
    Approximates the solution of a (mixed-)integer problem
    by iteratively computing a continuous(real) minimizer x followed by a
    rounding of the one variable x_i which is subsequently eliminated from the
    system, and the system is solved again ...
*/

class COMISODLLEXPORT MISolver
{
public:
   
  // typedefs
  typedef gmm::csc_matrix< double >       CSCMatrix;
  typedef std::vector< double >           Vecd;
  typedef std::vector< int >              Veci;
  typedef std::vector< unsigned int >     Vecui;

  // gmm Column and ColumnIterator of CSC matrix
  typedef gmm::linalg_traits< CSCMatrix >::const_sub_col_type Col;
  typedef gmm::linalg_traits< Col >::const_iterator           ColIter;


  /// default Constructor
  MISolver();

  /// Destructor
  ~MISolver() {}


  /// Compute greedy approximation to a mixed integer problem.
	/** @param _A symmetric positive semi-definite CSC matrix (Will be \b destroyed!)
	 *  @param _x vector holding solution at the end
   *  @param _rhs right hand side of system Ax=rhs (Will be \b destroyed!)
   *  @param _to_round vector with variable indices to round to integers
   *  @param _fixed_order specifies if _to_round indices shall be rounded in the
   *  given order (\b true) or be greedily selected (\b false)
	 *  */
  void solve(
    CSCMatrix& _A, 
    Vecd&      _x, 
    Vecd&      _rhs, 
    Veci&      _to_round,
    bool       _fixed_order = false );


  /// Compute greedy approximation to a mixed integer problem.
	/** @param _B mx(n+1) matrix with (still non-squared) equations of the energy,
   * including the right hand side (Will be \b destroyed!)
	 *  @param _x vector holding solution at the end
   *  @param _to_round vector with variable indices to round to integers
   *  @param _fixed_order specifies if _to_round indices shall be rounded in the
   *  given order (\b true) or be greedily selected (\b false)
	 *  */
  template<class CMatrixT>
  void solve( 
    CMatrixT& _B,
    Vecd&     _x,
    Veci&     _to_round,
    bool      _fixed_order = false );

  /// show Qt-Options-Dialog for setting algorithm parameters
  /** Requires a Qt Application running and COMISO_GUI to be defined */
  void show_options_dialog();

  /** @name Get/Set functions for algorithm parameters 
	 * Besides being used by the Qt-Dialog these can also be called explicitly
   * to set parameters of the algorithm. */
	/*@{*/
	/// Shall an initial full solution be computed?
  void set_inital_full( bool _b) {        initial_full_solution_=_b;}
	/// Will an initial full solution be computed?
  bool get_inital_full()         { return initial_full_solution_;}

	/// Shall a final full solution be computed?
  void set_final_full( bool _b) {        final_full_solution_=_b;}
	/// Will a final full solution be computed?
  bool get_final_full()         { return final_full_solution_;}

  /// Shall direct (or greedy) rounding be used?
  void set_direct_rounding( bool _b) {        direct_rounding_=_b;}
  /// Will direct rounding be used?
  bool get_direct_rounding()         { return direct_rounding_;}

  /// Shall no rounding be performed?
  void set_no_rounding( bool _b) {        no_rounding_=_b;}
  /// Will no rounding be performed?
  bool get_no_rounding()         { return no_rounding_;}

  /// Set number of maximum Gauss-Seidel iterations
  void         set_local_iters( unsigned int _i) { max_local_iters_ = _i;}
  /// Get number of maximum Gauss-Seidel iterations
  unsigned int get_local_iters()                 { return max_local_iters_;}

  /// Set error threshold for Gauss-Seidel solver
  void   set_local_error( double _d) { max_local_error_ = _d;}
  /// Get error threshold for Gauss-Seidel solver
  double get_local_error()           { return max_local_error_;}

  /// Set number of maximum Conjugate Gradient iterations 
  void         set_cg_iters( unsigned int _i) { max_cg_iters_ = _i;}
  /// Get number of maximum Conjugate Gradient iterations 
  unsigned int get_cg_iters()                 { return max_cg_iters_;}

  /// Set error threshold for Conjugate Gradient
  void   set_cg_error( double _d) { max_cg_error_ = _d;}
  /// Get error threshold for Conjugate Gradient
  double get_cg_error()           { return max_cg_error_;}

  /// Set error threshold before full solution is computed
  void   set_full_error( double _d) { max_full_error_ = _d;}
  /// Get error threshold before full solution is computed
  double get_full_error()           { return max_full_error_;}

  /// Set noise level of algorithm. 0 - quiet, 1 - more noise, 2 - even more, 100 - all noise
  void         set_noise( unsigned int _i) { noisy_ = _i;}
  /// Get noise level of algorithm
  unsigned int get_noise()                 { return noisy_;}

  /// Set output statistics of solver
  void set_stats( bool _stats) { stats_ = _stats; }
  /// Get output statistics of solver
  bool get_stats( )            { return stats_; }
	/*@}*/

private:

  /// Copy constructor (not used)
  MISolver(const MISolver& _rhs);

  /// Assignment operator (not used)
  MISolver& operator=(const MISolver& _rhs);

  // parameters used by the MiSo
  bool initial_full_solution_;
  bool final_full_solution_;

  bool direct_rounding_;
  bool no_rounding_;

  unsigned int max_local_iters_;
  double       max_local_error_;
  unsigned int max_cg_iters_;
  double       max_cg_error_;
  double       max_full_error_;
  unsigned int noisy_;
  bool         stats_;

  friend class ACG::MISolverDialog;
};


//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_MISOLVER_C)
#define ACG_MISOLVER_TEMPLATES
#include "MISolverT.cc"
#endif
//=============================================================================
#endif // ACG_MISOLVER_HH defined
//=============================================================================
