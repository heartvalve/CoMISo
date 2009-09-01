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



//=============================================================================
//
//  CLASS MISolver - IMPLEMENTATION
//
//=============================================================================

#define ACG_MISOLVER_C
//== INCLUDES =================================================================

#include "MISolver.hh"

//== NAMESPACES ===============================================================

namespace ACG {

//== IMPLEMENTATION ==========================================================
/*
template< class VecT>
void 
MISolver::solve( 
    gmm::col_matrix< VecT >& _B,
    VecT&                    _x,
    Veci&                    _to_round,
    bool                     _fixed_order = false )
{
  typedef typename gmm::row_matrix< VecT > RMatrixT;
  typedef typename gmm::col_matrix< VecT > CMatrixT;
  // some statistics
  int n_local = 0;
  int n_cg    = 0;
  int n_full  = 0;

  Veci to_round(_to_round);
  // if the order is not fixed, uniquify the indices
  if( !_fixed_order)
  {
    // copy to round vector and make it unique
    std::sort(to_round.begin(), to_round.end());
    Veci::iterator last_unique;
    last_unique = std::unique(to_round.begin(), to_round.end());
    int r = last_unique - to_round.begin();
    to_round.resize( r);
  }

  // setup quadratic csc system Ax=rhs
  CSCMatrix A;
  Vecd rhs;
  setup_quadratic_csc_system(_B, A, _x, rhs );

  // initalize old indices
  Vecui old_idx(rhs.size());
  for(unsigned int i=0; i<old_idx.size(); ++i)
    old_idx[i] = i;

  // Setup Cholmod solver used for full solution
  ACG::CholmodSolver chol;

  if( initial_full_solution_ || direct_rounding_)
  {
    if( noisy_ > 2) std::cerr << "initial full solution" << std::endl;
    chol.calc_system_gmm(A);
    chol.solve(_x, rhs);

    ++n_full;
  }

  if( no_rounding_ || _to_round.empty() )
    return;

  // preconditioner for CG
  gmm::identity_matrix PS, PR;
  // neighbors for local optimization
  Vecui neigh_i;

  // Vector for reduced solution
  Vecd xr(_x);

  // direct rounding?
  if( direct_rounding_ && _fixed_order )
    std::cerr << "Rounding order is fixed, direct rounding does not make sense!" << std::endl;
  if( direct_rounding_ && !_fixed_order )
  {
    int n_ints = to_round.size();
    int n_syst = _B.ncols();
    int rhs_idx= n_syst-1;

    // get rounded values and setup constraint matrix
    Veci c_elim; //columns to be eliminated
    RMatrixT C( n_ints, n_syst);
    n_ints = 0;
    for( uint j=0; j < to_round.size(); ++j)
    {
      if( to_round[j] > -1)
      {
        int cur_idx = to_round[j];
        double rnd_x = ROUND(xr[cur_idx]);
        C( j, cur_idx ) = -1;
        C( j, rhs_idx ) = rnd_x;
        c_elim.push_back( cur_idx);
        ++n_ints;
      }
    }
    if( to_round.size() != n_ints)
      gmm::resize( C, n_ints, n_syst);

    // eliminate all from _B (only economic when rounding very many)
    Veci empty_round_vector;
    Veci new_idx;
    eliminate_conditions( C, _B, empty_round_vector, c_elim, new_idx);

    // form new A
    setup_quadratic_csc_system( _B, A, _x, rhs);

    // perform final solve
    chol.calc_system_gmm(A);
    chol.solve(_x, rhs);
    ++n_full;
 
    // restore eliminated variables
    restore_eliminated_vars( C, _x, c_elim, new_idx);
   }
  if( !direct_rounding_)
   
  // loop until solution computed
  for(unsigned int i=0; i<to_round.size(); ++i)
  {
    if( noisy_ > 0)
    {
      std::cerr << "Integer DOF's left: " << to_round.size()-(i+1) << " ";
      if( noisy_ > 1)
        std::cerr << "residuum_norm: " << gmm::residuum_norm( A, xr, rhs) << std::endl;
    }

    // index to eliminate
    unsigned int i_best = 0;

    if( _fixed_order ) // if order is fixed, simply get next index from to_round
    {
      i_best = to_round[i];
    }
    else               // else search for best rounding candidate
    {
      // find index yielding smallest rounding error
      double       r_best = FLT_MAX;
      for(unsigned int j=0; j<to_round.size(); ++j)
      {
        if( to_round[j] != -1)
        {
          int cur_idx = to_round[j];
          double rnd_error = fabs( ROUND(xr[cur_idx]) - xr[cur_idx]);
          if( rnd_error < r_best)
          {
            i_best = cur_idx;
            r_best = rnd_error;
          }
        }
      }
    }

    // store rounded value
    double rnd_x = ROUND(xr[i_best]);
    _x[ old_idx[i_best] ] = rnd_x;

    // compute neighbors
    neigh_i.clear();
    Col col = gmm::mat_const_col(A, i_best);
    ColIter it  = gmm::vect_const_begin( col);
    ColIter ite = gmm::vect_const_end  ( col);
    for(; it!=ite; ++it)
      neigh_i.push_back(it.index());

    // eliminate var
    gmm::eliminate_var(i_best, rnd_x, A, xr, rhs);
    gmm::eliminate_var_idx(i_best, to_round);

    // update old_idx
    old_idx.erase( old_idx.begin()+i_best );

    if( direct_rounding_ && !_fixed_order) continue;

    // current error
    double cur_error = FLT_MAX;

    // compute new solution
    unsigned int n_its = max_local_iters_;
    if(max_local_iters_ > 0)
    {
      if( noisy_ > 2)std::cerr << "use local iteration ";

      n_its = gmm::gauss_seidel_local(A, xr, rhs, neigh_i, max_local_iters_, max_local_error_);
      ++n_local;
    }

    // update error bound
    // if gauss_seidel did not reach max iters then error must be less than
    // tolerance (max_local_error_)
    if( n_its < max_local_iters_) cur_error = max_local_error_;

    if( cur_error > max_cg_error_)
    {
      if( noisy_ > 2) std::cerr << ", cg ";

      gmm::iteration iter(max_cg_error_);
      iter.set_maxiter(max_cg_iters_);
      iter.set_noisy(std::max( int(noisy_)-4, int(0)));
      // conjugate gradient
      if( max_cg_iters_ > 0)
      {
        gmm::cg( A, xr, rhs, PS, PR, iter);
        cur_error = iter.get_res();
        ++n_cg;
      }
    }

    if( cur_error > max_full_error_ )
    {
      if( noisy_ > 2)std::cerr << ", full ";

      if( gmm::mat_ncols( A) > 0)
      {
        chol.calc_system_gmm(A);
        chol.solve(xr,rhs);

        ++n_full;
      }
    }

    if( noisy_ > 2)std::cerr << std::endl;
  }

  // final full solution?
  if( final_full_solution_ || direct_rounding_)
  {
    if( noisy_ > 2) std::cerr << "final full solution" << std::endl;

    if( gmm::mat_ncols( A) > 0)
    {
      chol.calc_system_gmm(A);
      chol.solve( xr, rhs);
      ++n_full;
    }
  }

  // store solution values to result vector
  for(unsigned int i=0; i<old_idx.size(); ++i)
  {
    _x[ old_idx[i] ] = xr[i];
  }

  // output statistics
  if( stats_)
  {
    std::cerr << "\t" << __FUNCTION__ << " *** Statistics of MiSo Solver ***\n";
    std::cerr << "\t\t Number of CG    iterations = " << n_cg << std::endl;
    std::cerr << "\t\t Number of LOCAL iterations = " << n_local << std::endl;
    std::cerr << "\t\t Number of FULL  iterations = " << n_full << std::endl;
    std::cerr << "\t\t Number of ROUNDING         = " << _to_round.size() << std::endl;
    std::cerr << std::endl;
  }

}

template< class CMatrixT >
void
MISolver::setup_quadratic_csc_system(
    CMatrixT&  _B
    CSCMatrix& _A,
    Vecd&      _x,
    Vecd&      _rhs )
{
  // clear _A
  _A.do_clear();

  // clear rhs
  _rhs.clear();

  unsigned int m = gmm::mat_nrows(_B);
  unsigned int n = gmm::mat_ncols(_B);

  // set up B transposed
  CMatrixT Bt( n ,m);
  gmm::copy( gmm::transposed( _B), Bt);

  // setup BtB
  CMatrixT BtB( n, n);
  gmm::mult( Bt, _B, BtB);

  // extract rhs
  _rhs.resize( n);
  gmm::copy( gmm::scaled(gmm::mat_const_col( BtB, n - 1),-1.0), rhs);
  _rhs.resize( n - 1);

  // resize BtB to only contain the actual system matrix (and not the rhs)
  gmm::resize( BtB, n - 1, n - 1);
  _x.resize( n - 1);

  // regularize if necessary
  //if(_reg_factor != 0.0)
  //  gmm::regularize_hack(BtB, _reg_factor);

  // BtB -> CSC
  _A.init_with_good_format( BtB);
}
*/
//-----------------------------------------------------------------------------



//=============================================================================
} // namespace ACG
//=============================================================================
