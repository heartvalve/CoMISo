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
//  CLASS Eigen_Tools - IMPLEMENTATION
//
//=============================================================================

#define COMISO_Eigen_TOOLS_C

//== INCLUDES =================================================================

#include "Eigen_Tools.hh"
#include <queue>
#include <CoMISo/Utils/StopWatch.hh>
#include <CoMISo/Utils/VSToolsT.hh>


//== NAMESPACES ===============================================================

namespace COMISO_EIGEN
{

//== IMPLEMENTATION ==========================================================


//-----------------------------------------------------------------------------

template<class MatrixT, class REALT, class INTT>
void get_ccs_symmetric_data( const MatrixT&      _mat,
                             const char          _uplo,
                             std::vector<REALT>& _values,
                             std::vector<INTT>&  _rowind,
                             std::vector<INTT>&  _colptr )
{
  // Assumes col major
  
   int m = _mat.innerSize();
   int n = _mat.outerSize();

   _values.resize( 0 );
   _rowind.resize( 0 );
   _colptr.resize( 0 );

   INTT iv( 0 );

   typedef typename MatrixT::InnerIterator It;

   switch ( _uplo )
   {
      case 'l':
      case 'L':
         // for all columns
         for ( int i=0; i<n; ++i )
         {
            _colptr.push_back( iv );

            // row it
            It it(_mat, i);

            for( ; it; ++it)
            {
              if( it.index() >= i )
              {
                _values.push_back( it.value());
                _rowind.push_back( it.index());
                ++iv;
              }
            }
         }
         _colptr.push_back( iv );
         break;

      case 'u':
      case 'U':
         // for all columns
         for ( int i=0; i<n; ++i )
         {
            _colptr.push_back( iv );

            // row it
            It it(_mat, i);

            for( ; it; ++it)
            {
              if( it.index() <= i )
              {
                _values.push_back( it.value());
                _rowind.push_back( it.index());
                ++iv;
              }
            }
         }
         _colptr.push_back( iv );
         break;

      case 'c':
      case 'C':
         // for all columns
         for ( int i=0; i<n; ++i )
         {
            _colptr.push_back( iv );

            // row it
            It it(_mat, i);

            for( ; it; ++it)
            {
              _values.push_back( it.value());
              _rowind.push_back( it.index());
              ++iv;
            }
         }
         _colptr.push_back( iv );
         break;

      default:
         std::cerr << "ERROR: parameter uplo must bei either 'U' or 'L' or 'C'!!!\n";
         break;
   }
}


//-----------------------------------------------------------------------------


// inspect the matrix: dimension, symmetry, zero_rows, zero_cols, nnz, max, min, max_abs, min_abs, NAN, INF
template<class MatrixT>
void inspect_matrix( const MatrixT& _A)
{


  std::cerr << "################### INSPECT MATRIX ##################\n";
  std::cerr << "#outer size  : " << _A.outerSize() << std::endl;
  std::cerr << "#inner size  : " << _A.innerSize() << std::endl;
  std::cerr << "#rows        : " << _A.rows() << std::endl;
  std::cerr << "#cols        : " << _A.cols() << std::endl;
  std::cerr << "#nonzeros    : " << _A.nonZeros() << std::endl;
  std::cerr << "#nonzeros/row: " << (double(_A.nonZeros())/double(_A.rows())) << std::endl;
  std::cerr << "symmetric    : " << is_symmetric( _A) << std::endl;

  MatrixT trans( _A.transpose());

  int zero_rows = 0;
  int zero_cols = 0;

  for(int i=0; i<_A.outerSize(); ++i)
  {
    typename MatrixT::InnerIterator it(_A, i); 
    if( !it) ++zero_rows;
  }

  for(int i=0; i<trans.outerSize(); ++i)
  {
    typename MatrixT::InnerIterator it(trans, i); 
    if( !it) ++zero_cols;
  }

  std::cerr << "zero rows    : " << zero_rows << std::endl;
  std::cerr << "zero cols    : " << zero_cols << std::endl;

  typedef typename MatrixT::Scalar Scalar;
  Scalar vmin     = std::numeric_limits<Scalar>::max();
  Scalar vmax     = std::numeric_limits<Scalar>::min();
  Scalar vmin_abs = std::numeric_limits<Scalar>::max();
  Scalar vmax_abs = 0;

  int n_nan = 0;
  int n_inf = 0;
  
  // inspect elements
  for(int i=0; i<_A.outerSize(); ++i)
  {
    typename MatrixT::InnerIterator it( _A, i);
    
    for(; it ; ++it)
    {
      if( it.value() < vmin ) vmin = it.value();
      if( it.value() > vmax ) vmax = it.value();

      if( fabs(it.value()) < vmin_abs) vmin_abs = fabs(it.value());
      if( fabs(it.value()) > vmax_abs) vmax_abs = fabs(it.value());

      if( std::isnan(it.value())) ++n_nan;
      if( std::isinf(it.value())) ++n_inf;
    }
  }
  
  std::cerr << "min  val     : " << vmin << std::endl;
  std::cerr << "max  val     : " << vmax << std::endl;
  std::cerr << "min |val|    : " << vmin_abs << std::endl;
  std::cerr << "max |val|    : " << vmax_abs << std::endl;
  std::cerr << "#nan         : " << n_nan << std::endl;
  std::cerr << "#inf         : " << n_inf << std::endl;
  
  std::cerr << "min eval     : " << "..." << std::endl;
  std::cerr << "max eval     : " << "..." << std::endl;
  std::cerr << "min|eval|    : " << "..." << std::endl;
  std::cerr << "max|eval|    : " << "..." << std::endl;
}

//-----------------------------------------------------------------------------


// symmetric ?
template<class MatrixT>
bool is_symmetric( const MatrixT& _A)
{
  typedef typename MatrixT::InnerIterator It;
  typedef typename MatrixT::Scalar Scalar;

  int nouter( _A.outerSize());
  int ninner( _A.innerSize());
  
  if( nouter != ninner )
    return false;

  bool symmetric(true);

  for( int c = 0; c < nouter; ++c)
  {
    for( It it(_A,c); it; ++it)
    {
      int r(it.index());

      Scalar val(it.value());

      // find diagonal partner element
      bool found(false);
      for( It dit(_A,r); dit; ++dit)
      {
        if( dit.index() < c )
        {}
        else if( dit.index() == c)
        {
          if( dit.value() == val) 
            found = true;
          break;
        }
        else 
        {
          break;
        }
      }
      if( !found) 
      {
        symmetric = false;
        break;
      }
    }
  }
  return symmetric;
}


//=============================================================================
} // namespace COMISO
//=============================================================================
