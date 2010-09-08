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


#include <CoMISo/Utils/StopWatch.hh>
#include <gmm/gmm.h>
#include <vector>
#include <CoMISo/Solver/ConstrainedSolver.hh>
#include <CoMISo/Solver/MISolver.hh>
#include <CoMISo/Solver/GMM_Tools.hh>


#define N  20
#define M  24
// previously 10, setback after TEST

template<class MatrixT, class CMatrixT>
void setup_system( MatrixT& _B, MatrixT& _Bt, CMatrixT& _A)
{
  gmm::resize(_B,  M, N);
  gmm::resize(_Bt, N, M);
  gmm::resize( _A, N, N);
  for( int i=0;i<M;++i)
    for( int j=0;j<N;++j)
      if( (rand()-0.8*RAND_MAX)/RAND_MAX> 0)
      _B(i,j) = round(((rand()-0.4*RAND_MAX)/RAND_MAX)*10.0);

  gmm::copy( gmm::transposed(_B), _Bt);
  gmm::mult( _Bt, _B, _A);
}

template<class Matrix, class Vector>
void get_rhs( Matrix& _A, Vector& _rhs)
{
  int size = gmm::mat_nrows(_A);
  _rhs.clear();
  _rhs.resize(size);
  gmm::copy( gmm::scaled(gmm::mat_const_col(_A, size-1),-1), _rhs);
  _rhs.resize( size-1);
}

template<class MatrixT>
void setup_constraints( MatrixT& _C, double _percentage, int _n = 1)
{
  gmm::resize( _C, _n, N);
  for( int i=0;i<_n;++i)
    for( int j=0;j<N;++j)
    {
      double randnum = (double(rand())/double(RAND_MAX));
      if ( randnum < _percentage)
        _C( i,j) = -1;
      else if( randnum > (1.0-_percentage))
        _C( i,j) = 1;
      else
        _C( i,j) = 0;
    }
}

int main(void)
{ 
  bool print(true);
  // setup matrix
  gmm::col_matrix< gmm::wsvector< double > > A;
  gmm::col_matrix< gmm::wsvector< double > > B;
  gmm::col_matrix< gmm::wsvector< double > > Bt;
  std::vector< double > rhs;
  setup_system(B, Bt, A);

  get_rhs( A, rhs);

  gmm::resize( A, N-1, N-1);

  // setup constraints
  gmm::row_matrix< gmm::wsvector< double > > C;
  setup_constraints( C, 0.2, 10);
/*  gmm::resize(C, 2,N);
  C(0,0) = 1;
  C(0,1) = 1;
  C(0,2) = 1;
  C(0,3) = 1;
  C(1,0) = 1;
  C(1,3) = 1;
   std::cerr << "C orig " << C << std::endl;*/


  // make conditions independent
  std::vector< int > idx_to_round, c_elim(gmm::mat_nrows(C));
  COMISO::ConstrainedSolver cs;
  cs.make_constraints_independent( C, idx_to_round, c_elim);
  //std::cerr << "C indep " << C << std::endl;
  //c_elim[0]=1;
  //std::cerr << "c_elim " << c_elim << std::endl;

  std::vector< double > x(N-1);
  std::vector< double > rhs_bak(rhs);

  // setup reindexing
  std::vector< int > new_idx(N, -1);
  for( unsigned int i = 0; i < N; ++i)
    new_idx[i] = i;
  std::cerr << "c_elim " << c_elim << std::endl;
 
  // test csc remove 
  gmm::csc_matrix< double > Mcsc( N-1, N-1);
  gmm::copy(A, Mcsc);
  std::vector<double> evals(c_elim.size());
  std::cerr << "c_elim " << c_elim << std::endl;
  std::cerr << "CSC BEFORE " << Mcsc << std::endl;
  COMISO_GMM::eliminate_csc_vars(c_elim, evals, Mcsc, x, rhs);
  std::cerr << "ELIMINATED CSC" << Mcsc << std::endl;

  //     T *pr;        // values.
//     IND_TYPE *ir; // row indices.
//     IND_TYPE *jc; // column repartition on pr and ir.
/*  gmm::fill_random(A);
  std::cerr << " A " << A << std::endl;
  gmm::row_matrix< gmm::wsvector< double > > Arow(gmm::mat_nrows(A), gmm::mat_ncols(A));
  gmm::copy(A, Arow);
  std::cerr << " Arow " << Arow << std::endl;

  std::cerr << "nc " << Mcsc.nc << " nr " << Mcsc.nr << " Mcsc " << Mcsc << std::endl;
  typedef unsigned int uint;

  uint elemcnt=0;
  for( uint c = 0; c < Mcsc.nc; ++c)
  {
    uint n_el = Mcsc.jc[c+1]-Mcsc.jc[c];
    for( uint i = 0; i < n_el; ++i)
    {
      std::cerr << "col " << c << " row " << Mcsc.ir[elemcnt] << " elem " << Mcsc.pr[elemcnt] << std::endl;
      ++elemcnt;
    }
  }
  for( unsigned int i = 0; i < Mcsc.nc; ++i)
  {
    std::cerr << "pr " << Mcsc.pr[i] << " ir " << Mcsc.ir[i] << " jc " << Mcsc.jc[i] <<  std::endl;
  }
    std::cerr << "jc :";*/
/*  for( unsigned int i = 0; i < Mcsc.nc*Mcsc.nr; ++i)
    std::cerr << Mcsc.jc[i] << " ";
  std::cerr << std::endl;
*/

  







/*





  // test quadratic elimination
  //std::cerr << " A " << A << std::endl;
  COMISO::StopWatch sw;
  sw.start();
  gmm::eliminate_constraints(C, A, x, rhs, idx_to_round, c_elim, new_idx);

  std::cerr << "\t TIME QUADRATIC " << sw.stop()/1000.0 << std::endl;
  gmm::col_matrix< gmm::wsvector< double > > Aquad(A);
  //std::cerr << " rhs after elim " << rhs << std::endl;

  // WARNING TODO somehow memory leak when doing resize in eliminiate_vars... do it here
  // TODO also adapt "real" eliminate vars function in PHYSIM!!!! DO NOW!!
  //gmm::resize(Aquad, rhs.size(), rhs.size());
  std::cerr << " A after elim " << Aquad << std::endl;
  
  gmm::col_matrix< gmm::rsvector< double > > Bcol;

  //std::cerr << "c_elim " << c_elim << std::endl;
  //std::cerr << " C " << C << std::endl;
  gmm::row_matrix< gmm::wsvector<double> > Brow(gmm::mat_nrows(B), gmm::mat_ncols(B));
  gmm::copy(B,Brow);
  COMISO::StopWatch sw1;
  COMISO::StopWatch sw2;
  COMISO::StopWatch sw3;
  double time1=0.0, time2=0.0, time3=0.0;
  gmm::col_matrix< gmm::wsvector<double > > Btt;
  sw.start();

  sw1.start();
  cs.eliminate_conditions( C,  Brow, idx_to_round, c_elim, new_idx, Bcol);
  time1 = sw1.stop();
  sw2.start();
  gmm::resize(Btt, Bcol.ncols(), Bcol.nrows());
  gmm::copy(gmm::transposed(Bcol), Btt);
  gmm::resize(A, Bcol.ncols(), Bcol.ncols());
  time2=sw2.stop();
  sw3.start();
  gmm::mult(Btt,Bcol,A);
  time3=sw3.stop();
  std::cerr << "\t TIME FACTORED " << sw.stop()/1000.0 << std::endl;
  std::cerr << "\t eliminiate took " << time1/1000.0 << std::endl;
  std::cerr << "\t copy took       " << time2/1000.0 << std::endl;
  std::cerr << "\t muliply took    " << time3/1000.0 << std::endl;
  gmm::resize(A, gmm::mat_nrows(A)-1, gmm::mat_ncols(A)-1);
  std::cerr << " A is " << gmm::mat_nrows(A) << " x " << gmm::mat_ncols(A) << " Aquad is " << gmm::mat_nrows(Aquad) << " x " << gmm::mat_ncols(Aquad) << std::endl;
  gmm::add(gmm::scaled(A,-1),Aquad);
  std::cerr << " A " << A << std::endl;
  std::cerr << " Aquad " << Aquad << std::endl;
  std::cerr << "MATRIX NORM " << gmm::mat_norminf(Aquad) << std::endl;










  */





/*  std::cerr << "Standard elimination with deletion between every..." << std::endl;
  std::cerr << " A " << A << std::endl;
  for( int i=0;i<nCrows-1;++i)
  {
    int idx=c_elim[i];
    if( idx != -1)
    {
      gmm::wsvector<double> constraint(C.row(i));
      typedef gmm::linalg_traits<gmm::wsvector<double> >::iterator CIter;
      CIter it     = gmm::vect_begin( constraint );
      CIter it_end = gmm::vect_end( constraint );
      double val = constraint[idx];

      for ( ; it != it_end; ++it )
        (*it) /= val;
      constraint[idx]=0;
      std::cerr << "idx " << idx << " constraint " << constraint << " rhs " << Crhs[i] << std::endl;
      gmm::eliminate_condition( idx, constraint, Crhs[i], A, x, rhs);
    }
  }
  std::cerr << " A " << A << std::endl;

  std::cerr << "It took " << sw.stop()/1000.0 << std::endl;

  gmm::col_matrix< gmm::rsvector< double > > Bcol;
  std::vector<int>                          new_idxx;
  std::vector<int> itoround;
  cs.eliminate_conditions( C,  B, itoround, c_elim, new_idxx, Bcol);
  gmm::resize(B, Bcol.nrows(), Bcol.ncols());
  gmm::copy(Bcol,B);
  setup_system(B, Bt, A);
  std::cerr << " A " << A << std::endl;
*/

/*  // setup constraints

  eliminate_condition( 1, gmm::mat_const_row(C,0), 1, A, x, rhs);

  std::cerr << "Standard elimination withOUT deletion between every..." << std::endl;
  std::cerr << A << std::endl;

  std::cerr << rhs_bak << std::endl << rhs << std::endl;
*/
  return -1;
}

