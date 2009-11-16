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

#include <gmm/gmm.h>
#include <float.h>

#define ROUND(x) ((x)<0?int((x)-0.5):int((x)+0.5))

namespace COMISO {



// Constructor
MISolver::MISolver() 
{
  // default parameters
  initial_full_solution_ = true;
  final_full_solution_   = true;

  direct_rounding_ = false;
  no_rounding_     = false;

  max_local_iters_ = 10000;
  max_local_error_ = 1e-3;
  max_cg_iters_    = 20;
  max_cg_error_    = 1e-3;
  max_full_error_   = 1e-1;
  noisy_ = 0;
  stats_ = true;
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

  // initalize old indices
  Veci old_idx(_rhs.size());
  for(unsigned int i=0; i<old_idx.size(); ++i)
    old_idx[i] = i;

  // Setup Cholmod solver used for full solution
  COMISO::CholmodSolver chol;

  if( initial_full_solution_ || direct_rounding_)
  {
    if( noisy_ > 2) std::cerr << "initial full solution" << std::endl;
    chol.calc_system_gmm(_A);
    chol.solve(_x, _rhs);

    ++n_full;
  }

 if( no_rounding_)
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

 // direct rounding?
  if( direct_rounding_)
  {
    Vecui elim_i;
    Vecd  elim_v;
    for( unsigned int i=0; i < to_round.size(); ++i)
    {
      _x[to_round[i]] = ROUND(xr[to_round[i]]);
      elim_i.push_back(to_round[i]);
      elim_v.push_back(_x[to_round[i]]);
      // update old idx
      old_idx[to_round[i]] = -1;
    }
    Veci::iterator new_end = std::remove( old_idx.begin(), old_idx.end(), -1);
    old_idx.resize( new_end-old_idx.begin());
    // eliminate vars from linear system
    gmm::eliminate_csc_vars2( elim_i, elim_v, _A, xr, _rhs);
  }

  else

  // loop until solution computed
  for(unsigned int i=0; i<to_round.size(); ++i)
  {
    if( noisy_ > 0)
    {
      std::cerr << "Integer DOF's left: " << to_round.size()-(i+1) << " ";
      if( noisy_ > 1)
        std::cerr << "residuum_norm: " << gmm::residuum_norm( _A, xr, _rhs) << std::endl;
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
    Col col = gmm::mat_const_col(_A, i_best);
    ColIter it  = gmm::vect_const_begin( col);
    ColIter ite = gmm::vect_const_end  ( col);
    for(; it!=ite; ++it)
      if(it.index() < i_best)
        neigh_i.push_back(it.index());
      else
        if(it.index() > i_best)
          neigh_i.push_back(it.index()-1);


    // eliminate var
    gmm::eliminate_var(i_best, rnd_x, _A, xr, _rhs);
    gmm::eliminate_var_idx(i_best, to_round);

    // update old_idx
    old_idx.erase( old_idx.begin()+i_best );

    //if( direct_rounding_ && !_fixed_order) continue;

    // current error
    double cur_error = FLT_MAX;

    // compute new solution
    unsigned int n_its = max_local_iters_;
    if(max_local_iters_ > 0)
    {
      if( noisy_ > 2)std::cerr << "use local iteration ";

      n_its = gmm::gauss_seidel_local(_A, xr, _rhs, neigh_i, max_local_iters_, max_local_error_);
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
        gmm::cg( _A, xr, _rhs, PS, PR, iter);
        cur_error = iter.get_res();
        ++n_cg;
      }
    }

    if( cur_error > max_full_error_ )
    {
      if( noisy_ > 2)std::cerr << ", full ";

      if( gmm::mat_ncols( _A) > 0)
      {
        chol.calc_system_gmm(_A);
        chol.solve(xr,_rhs);

        ++n_full;
      }
    }

    if( noisy_ > 2)std::cerr << std::endl;
  }

  // final full solution?
  if( _to_round.size() != 0)
  if( final_full_solution_ || direct_rounding_)
  {
    if( noisy_ > 2) std::cerr << "final full solution" << std::endl;

    if( gmm::mat_ncols( _A) > 0)
    {
      chol.calc_system_gmm(_A);
      chol.solve( xr, _rhs);
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
