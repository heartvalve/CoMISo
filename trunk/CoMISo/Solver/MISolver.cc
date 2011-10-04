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



#include "MISolver.hh"

#ifdef QT4_FOUND
#include <CoMISo/QtWidgets/MISolverDialogUI.hh>
#endif

#include <CoMISo/Utils/StopWatch.hh>

#include <gmm/gmm.h>
#include <float.h>

// hack for testing only
#include "SparseQRSolver.hh"
#include "UMFPACKSolver.hh"

#define ROUND(x) ((x)<0?int((x)-0.5):int((x)+0.5))

namespace COMISO {



// Constructor
MISolver::MISolver() 
{
  // default parameters
  initial_full_solution_ = true;
  iter_full_solution_    = true;
  final_full_solution_   = true;

  direct_rounding_   = false;
  no_rounding_       = false;
  multiple_rounding_ = true;

  max_local_iters_ = 100000;
  max_local_error_ = 1e-3;
  max_cg_iters_    = 50;
  max_cg_error_    = 1e-3;

  multiple_rounding_threshold_ = 0.5;

  noisy_ = 0;
  stats_ = true;

  use_constraint_reordering_ = true;
}


//-----------------------------------------------------------------------------


void 
MISolver::solve( 
    CSCMatrix& _A, 
    Vecd&      _x, 
    Vecd&      _rhs, 
    Veci&      _to_round,
    bool       _fixed_order )
{
  // nothing to solve?
  if( gmm::mat_ncols(_A) == 0 || gmm::mat_nrows(_A) == 0)
    return;

  if( no_rounding_ || _to_round.size() == 0)
    solve_no_rounding( _A, _x, _rhs);
  else
    if( direct_rounding_)
      solve_direct_rounding( _A, _x, _rhs, _to_round);
    else
      if( multiple_rounding_)
	solve_multiple_rounding( _A, _x, _rhs, _to_round);
      else
	solve_iterative( _A, _x, _rhs, _to_round, _fixed_order);
}


//-----------------------------------------------------------------------------


void 
MISolver::solve_no_rounding( 
    CSCMatrix& _A, 
    Vecd&      _x, 
    Vecd&      _rhs )
{
  chol_.calc_system_gmm(_A);
  chol_.solve(_x, _rhs);
}


//-----------------------------------------------------------------------------


void 
MISolver::solve_direct_rounding( 
    CSCMatrix& _A, 
    Vecd&      _x, 
    Vecd&      _rhs, 
    Veci&      _to_round)
{
  Veci to_round(_to_round);
  // copy to round vector and make it unique
  std::sort(to_round.begin(), to_round.end());
  Veci::iterator last_unique;
  last_unique = std::unique(to_round.begin(), to_round.end());
  int r = last_unique - to_round.begin();
  to_round.resize( r);

  // initalize old indices
  Veci old_idx(_rhs.size());
  for(unsigned int i=0; i<old_idx.size(); ++i)
    old_idx[i] = i;
  chol_.calc_system_gmm(_A);
  chol_.solve(_x, _rhs);

  // check solver performance (only for testing!!!)
  {
    StopWatch sw;

    // performance comparison code
    {
      sw.start();
      COMISO::SparseQRSolver spqr;
      spqr.calc_system_gmm(_A);
      std::cerr << "SparseQR factor took: " << sw.stop()/1000.0 << "s\n";
      Vecd x2(_x);
      sw.start();
      spqr.solve(x2,_rhs);
      std::cerr << "SparseQR solve took: " << sw.stop()/1000.0 << "s\n";
      Vecd res(_x);
      gmm::add(_x,gmm::scaled(x2,-1.0),res);
      std::cerr << "DIFFERENCE IN RESULT: " << gmm::vect_norm2(res) << std::endl;
    }

    // performance comparison code
    {
      sw.start();
      COMISO::UMFPACKSolver umf;
      umf.calc_system_gmm(_A);
      std::cerr << "UMFPack factor took: " << sw.stop()/1000.0 << "s\n";
      Vecd x3(_x);
      sw.start();
      umf.solve(x3,_rhs);
      std::cerr << "UMFPack solve took: " << sw.stop()/1000.0 << "s\n";
      Vecd res2(_x);
      gmm::add(_x,gmm::scaled(x3,-1.0),res2);
      std::cerr << "UMFPACK DIFFERENCE IN RESULT: " << gmm::vect_norm2(res2) << std::endl;
    }

    // performance comparison code
    {
      sw.start();
      COMISO::CholmodSolver chol;
      chol.calc_system_gmm(_A);
      std::cerr << "Choldmod factor took: " << sw.stop()/1000.0 << "s\n";
      Vecd x4(_x);
      sw.start();
      chol.solve(x4,_rhs);
      std::cerr << "Choldmod solve took: " << sw.stop()/1000.0 << "s\n";
      Vecd res(_x);
      gmm::add(_x,gmm::scaled(x4,-1.0),res);
      std::cerr << "DIFFERENCE IN RESULT: " << gmm::vect_norm2(res) << std::endl;
    }
  }

  // round and eliminate variables
  Vecui elim_i;
  Vecd  elim_v;
  for( unsigned int i=0; i < to_round.size(); ++i)
  {
    _x[to_round[i]] = ROUND(_x[to_round[i]]);
    elim_i.push_back(to_round[i]);
    elim_v.push_back(_x[to_round[i]]);
    // update old idx
    old_idx[to_round[i]] = -1;
  }

  Veci::iterator new_end = std::remove( old_idx.begin(), old_idx.end(), -1);
  old_idx.resize( new_end-old_idx.begin());
  // eliminate vars from linear system
  Vecd xr(_x);
  COMISO_GMM::eliminate_csc_vars2( elim_i, elim_v, _A, xr, _rhs);

  // std::cerr << "size A: " << gmm::mat_nrows(_A) << " " << gmm::mat_ncols(_A) 
  // 	    << std::endl;
  // std::cerr << "size x  : " << xr.size() << std::endl;
  // std::cerr << "size rhs: " << _rhs.size() << std::endl;

  // final full solution
  if( gmm::mat_ncols( _A) > 0)
  {
    //    chol_.update_system_gmm(_A);
    chol_.calc_system_gmm(_A);
    chol_.solve( xr, _rhs);
  }

  // store solution values to result vector
  for(unsigned int i=0; i<old_idx.size(); ++i)
  {
    _x[ old_idx[i] ] = xr[i];
  }
}


//-----------------------------------------------------------------------------


void 
MISolver::solve_iterative( 
    CSCMatrix& _A, 
    Vecd&      _x, 
    Vecd&      _rhs, 
    Veci&      _to_round,
    bool       _fixed_order )
{
  // StopWatch
  COMISO::StopWatch sw;
  double time_search_next_integer = 0;

  // some statistics
  n_local_ = 0;
  n_cg_    = 0;
  n_full_  = 0;

  // reset cholmod step flag
  cholmod_step_done_ = false;

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

  // initalize old indices
  Veci old_idx(_rhs.size());
  for(unsigned int i=0; i<old_idx.size(); ++i)
    old_idx[i] = i;

  if( initial_full_solution_)
  {
    if( noisy_ > 2) std::cerr << "initial full solution" << std::endl;
    chol_.calc_system_gmm(_A);
    chol_.solve(_x, _rhs);

    cholmod_step_done_ = true;

    ++n_full_;
  }

  // neighbors for local optimization
  Vecui neigh_i;

  // Vector for reduced solution
  Vecd xr(_x);

  // loop until solution computed
  for(unsigned int i=0; i<to_round.size(); ++i)
  {
    if( noisy_ > 0)
    {
      std::cerr << "Integer DOF's left: " << to_round.size()-(i+1) << " ";
      if( noisy_ > 1)
        std::cerr << "residuum_norm: " << COMISO_GMM::residuum_norm( _A, xr, _rhs) << std::endl;
    }

    // index to eliminate
    unsigned int i_best = 0;

    // position in round vector
    unsigned int tr_best = 0;

    if( _fixed_order ) // if order is fixed, simply get next index from to_round
    {
      i_best = to_round[i];
    }
    else               // else search for best rounding candidate
    {
      sw.start();
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
            i_best  = cur_idx;
            r_best  = rnd_error;
	    tr_best = j;
          }
        }
      }
      time_search_next_integer += sw.stop();
    }

    // store rounded value
    double rnd_x = ROUND(xr[i_best]);
    _x[ old_idx[i_best] ] = rnd_x;

    // compute neighbors
    neigh_i.clear();
    Col col = gmm::mat_const_col(_A, i_best);
    ColIter it  = gmm::vect_const_begin( col);
    ColIter ite = gmm::vect_const_end  ( col);
    for(; it!=ite; ++it)
      if(it.index() != i_best)
        neigh_i.push_back(it.index());

    // eliminate var
    COMISO_GMM::fix_var_csc_symmetric(i_best, rnd_x, _A, xr, _rhs);
    to_round[tr_best] = -1;

    // 3-stage update of solution w.r.t. roundings
    // local GS / CG / SparseCholesky
    update_solution( _A, xr, _rhs, neigh_i);
  }

  // final full solution?
  if( final_full_solution_)
  {
    if( noisy_ > 2) std::cerr << "final full solution" << std::endl;

    if( gmm::mat_ncols( _A) > 0)
    {
      if(cholmod_step_done_)
	chol_.update_system_gmm(_A);
      else
	chol_.calc_system_gmm(_A);

      chol_.solve( xr, _rhs);
      ++n_full_;
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
    std::cerr << "\t\t Number of CG    iterations  = " << n_cg_ << std::endl;
    std::cerr << "\t\t Number of LOCAL iterations  = " << n_local_ << std::endl;
    std::cerr << "\t\t Number of FULL  iterations  = " << n_full_ << std::endl;
    std::cerr << "\t\t Number of ROUNDING          = " << _to_round.size() << std::endl;
    std::cerr << "\t\t time searching next integer = " << time_search_next_integer / 1000.0 <<"s\n";
    std::cerr << std::endl;
  }
}



//-----------------------------------------------------------------------------


void 
MISolver::update_solution( 
    CSCMatrix& _A, 
    Vecd&      _x, 
    Vecd&      _rhs, 
    Vecui&     _neigh_i )
{
  // set to not converged
  bool converged = false;

  // compute new solution
  if(max_local_iters_ > 0)
  {
    if( noisy_ > 2)std::cerr << "use local iteration ";

    int    n_its     = max_local_iters_;
    double tolerance = max_local_error_;
    converged = siter_.gauss_seidel_local(_A, _x, _rhs, _neigh_i, n_its, tolerance);

    ++n_local_;
  }


  // conjugate gradient
  if( !converged && max_cg_iters_ > 0)
  {
    if( noisy_ > 2) std::cerr << ", cg ";

    int max_cg_iters = max_cg_iters_;
    double tolerance = max_cg_error_;
    converged = siter_.conjugate_gradient(_A, _x,_rhs, max_cg_iters, tolerance);

    if( noisy_ > 3) 
      std::cerr << "( converged " << converged << " "
		<< " iters " << max_cg_iters   << " "
		<< " res_norm " << tolerance << std::endl;
    ++n_cg_;
  }

  if(!converged && iter_full_solution_)
  {
    if( noisy_ > 2)std::cerr << ", full ";

    if( gmm::mat_ncols( _A) > 0)
    {
      if(cholmod_step_done_)
	chol_.update_system_gmm(_A);
      else
      {
	chol_.calc_system_gmm(_A);
	cholmod_step_done_ = true;
      }
      chol_.solve(_x,_rhs);

      ++n_full_;
    }
  }

  if( noisy_ > 2)std::cerr << std::endl;
}

//-----------------------------------------------------------------------------


void 
MISolver::solve_multiple_rounding( 
    CSCMatrix& _A, 
    Vecd&      _x, 
    Vecd&      _rhs, 
    Veci&      _to_round )
{
  // StopWatch
  COMISO::StopWatch sw;
  double time_search_next_integer = 0;

  // some statistics
  n_local_ = 0;
  n_cg_    = 0;
  n_full_  = 0;

  // reset cholmod step flag
  cholmod_step_done_ = false;

  Veci to_round(_to_round);
  // copy to round vector and make it unique
  std::sort(to_round.begin(), to_round.end());
  Veci::iterator last_unique;
  last_unique = std::unique(to_round.begin(), to_round.end());
  int r = last_unique - to_round.begin();
  to_round.resize( r);

  // initalize old indices
  Veci old_idx(_rhs.size());
  for(unsigned int i=0; i<old_idx.size(); ++i)
    old_idx[i] = i;

  if( initial_full_solution_)
  {
    if( noisy_ > 2) std::cerr << "initial full solution" << std::endl;
    chol_.calc_system_gmm(_A);
    chol_.solve(_x, _rhs);

    cholmod_step_done_ = true;

    ++n_full_;
  }

  // neighbors for local optimization
  Vecui neigh_i;

  // Vector for reduced solution
  Vecd xr(_x);

  // loop until solution computed
  for(unsigned int i=0; i<to_round.size(); ++i)
  {
    if( noisy_ > 0)
    {
      std::cerr << "Integer DOF's left: " << to_round.size()-(i+1) << " ";
      if( noisy_ > 1)
        std::cerr << "residuum_norm: " << COMISO_GMM::residuum_norm( _A, xr, _rhs) << std::endl;
    }

    // position in round vector
    std::vector<int> tr_best;

    sw.start();

    RoundingSet rset;
    rset.set_threshold(multiple_rounding_threshold_);

    // find index yielding smallest rounding error
    for(unsigned int j=0; j<to_round.size(); ++j)
    {
      if( to_round[j] != -1)
      {
	int cur_idx = to_round[j];
	double rnd_error = fabs( ROUND(xr[cur_idx]) - xr[cur_idx]);

	rset.add(j, rnd_error);
      }
    }

    rset.get_ids( tr_best);

    time_search_next_integer += sw.stop();
  
    // nothing more to do?
    if( tr_best.size() == 0)
      break;

    if( noisy_ > 5)
      std::cerr << "round " << tr_best.size() << " variables simultaneously\n";

    // clear neigh for local update
    neigh_i.clear();

    for(unsigned int j = 0; j<tr_best.size(); ++j)
    {
      int i_cur = to_round[tr_best[j]];

      // store rounded value
      double rnd_x = ROUND(xr[i_cur]);
      _x[ old_idx[i_cur] ] = rnd_x;

      // compute neighbors
      Col col = gmm::mat_const_col(_A, i_cur);
      ColIter it  = gmm::vect_const_begin( col);
      ColIter ite = gmm::vect_const_end  ( col);
      for(; it!=ite; ++it)
	if(it.index() != (unsigned int)i_cur)
	  neigh_i.push_back(it.index());

      // eliminate var
      COMISO_GMM::fix_var_csc_symmetric( i_cur, rnd_x, _A, xr, _rhs);
      to_round[tr_best[j]] = -1;
    }

    // 3-stage update of solution w.r.t. roundings
    // local GS / CG / SparseCholesky
    update_solution( _A, xr, _rhs, neigh_i);
  }

  // final full solution?
  if( final_full_solution_)
  {
    if( noisy_ > 2) std::cerr << "final full solution" << std::endl;

    if( gmm::mat_ncols( _A) > 0)
    {
      if(cholmod_step_done_)
	chol_.update_system_gmm(_A);
      else
	chol_.calc_system_gmm(_A);

      chol_.solve( xr, _rhs);
      ++n_full_;
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
    std::cerr << "\t\t Number of CG    iterations  = " << n_cg_ << std::endl;
    std::cerr << "\t\t Number of LOCAL iterations  = " << n_local_ << std::endl;
    std::cerr << "\t\t Number of FULL  iterations  = " << n_full_ << std::endl;
    std::cerr << "\t\t Number of ROUNDING          = " << _to_round.size() << std::endl;
    std::cerr << "\t\t time searching next integer = " << time_search_next_integer / 1000.0 <<"s\n";
    std::cerr << std::endl;
  }
}


//----------------------------------------------------------------------------


void 
MISolver::
show_options_dialog()
{
#ifdef QT4_FOUND
  MISolverDialog* pd = new MISolverDialog(*this);
  pd->exec();
#endif
}


// end namespace COMISO
}// ----------------------