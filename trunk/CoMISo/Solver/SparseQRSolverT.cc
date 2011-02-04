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



#define COMISO_SPARSE_QR_SOLVER_TEMPLATES_C

#include "SparseQRSolver.hh"


namespace COMISO {


template< class GMM_MatrixT>
bool SparseQRSolver::calc_system_gmm( const GMM_MatrixT& _mat)
{
//   std::vector<int>    colptr;
//   std::vector<int>    rowind;
//   std::vector<double> values;
    

    if(show_timings_) sw_.start();

    COMISO_GMM::get_ccs_symmetric_data( _mat,
					'c',
					values_, 
					rowind_, 
					colptr_ );
    
    if(show_timings_)
    {
      std::cerr << "Cholmod Timing GMM convert: " << sw_.stop()/1000.0 << "s\n";
    }

    return calc_system( colptr_, rowind_, values_);
}
  

//-----------------------------------------------------------------------------


template< class GMM_MatrixT>
bool SparseQRSolver::update_system_gmm( const GMM_MatrixT& _mat)
{
//   std::vector<int>    colptr;
//   std::vector<int>    rowind;
//   std::vector<double> values;
    
  COMISO_GMM::get_ccs_symmetric_data( _mat,
				      'c',
				      values_, 
				      rowind_, 
				      colptr_ );

    return update_system( colptr_, rowind_, values_);
}


//-----------------------------------------------------------------------------



template< class GMM_MatrixT, class GMM_MatrixT2, class GMM_MatrixT3, class IntT>
int
SparseQRSolver::
factorize_system_gmm( const GMM_MatrixT& _A, GMM_MatrixT2& _Q, GMM_MatrixT3& _R, std::vector<IntT>& _P)
{
  std::cerr << "factorize_system_gmm" << std::endl;
  // get dimensions
  int m = gmm::mat_nrows(_A);
  int n = gmm::mat_ncols(_A);

  // 1. _A -> cholmod_sparse A
  cholmod_sparse* AC(0);
  COMISO_GMM::gmm_to_cholmod(_A, AC, mp_cholmodCommon, 0, true);
  std::cerr << "gmm_to_cholmod finished" << std::endl;
  cholmod_print_sparse(AC, "AC", mp_cholmodCommon);

  // 2. factorize A -> Q,R,P
  UF_long econ = m;
  cholmod_sparse *Q, *R;
//  UF_long *P = new UF_long[n];
  UF_long *P;
  double rank = SuiteSparseQR<double>(ordering_, tolerance_, econ, AC, &Q, &R, &P, mp_cholmodCommon);
  std::cerr << "factorization finished" << std::endl;
  std::cerr << "rank: " << rank << std::endl;
  cholmod_print_sparse(Q, "Q", mp_cholmodCommon);

  // 3. convert Q,R,P -> _Q, _R, _P
  COMISO_GMM::cholmod_to_gmm(*Q, _Q);
  COMISO_GMM::cholmod_to_gmm(*R, _R);
  std::cerr << "cholmod_to_gmm finished" << std::endl;

  _P.clear(); _P.resize(n);
  for( int i=0; i<n; ++i)
    _P[i] = P[i];
  std::cerr << "coy vector finished" << std::endl;

  cholmod_l_free_sparse(&Q, mp_cholmodCommon);
  cholmod_l_free_sparse(&R, mp_cholmodCommon);
  std::cerr << "free1 finished" << std::endl;

  // TODO: alloc or free P ???
  cholmod_free(n, sizeof(UF_long), P, mp_cholmodCommon);
  std::cerr << "free2 finished" << std::endl;


  //// [Q,R,E] = qr(A), returning Q as a sparse matrix
//template <typename Entry> UF_long SuiteSparseQR     // returns rank(A) estimate
//(
//    int ordering,           // all, except 3:given treated as 0:fixed
//    double tol,
//    UF_long econ,
//    cholmod_sparse *A,      // m-by-n sparse matrix
//    // outputs
//    cholmod_sparse **Q,     // m-by-e sparse matrix where e=max(econ,rank(A))
//    cholmod_sparse **R,     // e-by-n sparse matrix
//    UF_long **E,            // permutation of 0:n-1, NULL if identity
//    cholmod_common *cc      // workspace and parameters
//) ;

  std::cerr << " ############## QR Factorization Info #############\n";
  std::cerr << " m: " << m << ", n: " << n << ", rank: " << rank << std::endl;

  return rank;
}

}
