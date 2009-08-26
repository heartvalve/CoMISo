
#include "cholmod.h"

#include "CholmodSolver.hh"


namespace ACG {

CholmodSolver::CholmodSolver()
{
    mp_cholmodCommon = new cholmod_common;
    cholmod_start( mp_cholmodCommon );

    mp_L = 0;
}


  //-----------------------------------------------------------------------------
  

CholmodSolver::~CholmodSolver()
{
    if( mp_L )
    {
	cholmod_free_factor( &mp_L, mp_cholmodCommon );
    }

    cholmod_finish( mp_cholmodCommon );
}
  

//-----------------------------------------------------------------------------


bool CholmodSolver::calc_system( const std::vector<int>&    _colptr, 
				 const std::vector<int>&    _rowind, 
				 const std::vector<double>& _values)
{
    colptr_ = _colptr;
    rowind_ = _rowind;
    values_ = _values;

    int n   = colptr_.size()-1;

    cholmod_sparse matA;

    matA.nrow = n;
    matA.ncol = n;
    matA.nzmax = _values.size();

    matA.p = &colptr_[0];
    matA.i = &rowind_[0];
    matA.x = &values_[0];
    matA.nz = 0;
    matA.z = 0;
    
    matA.stype = -1;
    matA.itype = CHOLMOD_INT;
    matA.xtype = CHOLMOD_REAL;
    matA.dtype = CHOLMOD_DOUBLE;
    matA.sorted = 1;
    matA.packed = 1;

    // clean up
    if( mp_L )
    {
	cholmod_free_factor( &mp_L, mp_cholmodCommon );
	mp_L = 0;
    }


    if( !(mp_L = cholmod_analyze( &matA, mp_cholmodCommon )) )
    {
	std::cout << "cholmod_analyze failed" << std::endl;
	return false;
    }

    if( !cholmod_factorize( &matA, mp_L, mp_cholmodCommon ) )
    {
	std::cout << "cholmod_factorize failed" << std::endl;
	return false;
    }


    return true;
}


  //-----------------------------------------------------------------------------

    
bool CholmodSolver::update_system( const std::vector<int>& _colptr, 
				   const std::vector<int>& _rowind, 
				   const std::vector<double>& _values )
{
    if( !mp_L )
	return false;

    colptr_ = _colptr;
    rowind_ = _rowind;
    values_ = _values;
    int n   = colptr_.size()-1;

    cholmod_sparse matA;

    matA.nrow = n;
    matA.ncol = n;
    matA.nzmax = _values.size();

    matA.p = &colptr_[0];
    matA.i = &rowind_[0];
    matA.x = &values_[0];
    matA.nz = 0;
    matA.z = 0;
    
    matA.stype = -1;
    matA.itype = CHOLMOD_INT;
    matA.xtype = CHOLMOD_REAL;
    matA.dtype = CHOLMOD_DOUBLE;
    matA.sorted = 1;
    matA.packed = 1;


    if( !cholmod_factorize( &matA, mp_L, mp_cholmodCommon ) )
    {
	std::cout << "cholmod_factorize failed" << std::endl;
	return false;
    }

    return true;
}


//-----------------------------------------------------------------------------
  

bool CholmodSolver::solve( double * _x, double * _b)
{
    const unsigned int n = colptr_.size() - 1;

    cholmod_dense *x, b;


    b.nrow = n;
    b.ncol = 1;
    b.nzmax = n;
    b.d = b.nrow;
    b.x = _b;
    b.z = 0;
    b.xtype = CHOLMOD_REAL;
    b.dtype = CHOLMOD_DOUBLE;

    if( !(x = cholmod_solve( CHOLMOD_A, mp_L, &b, mp_cholmodCommon )) )
    {
	std::cout << "cholmod_solve failed" << std::endl;
	return false;
    }
    
    for( unsigned int i = 0; i < n; ++i )
	_x[i] = ((double*)x->x)[i];

    cholmod_free_dense( &x, mp_cholmodCommon );

    return true;
}


//-----------------------------------------------------------------------------


  
}
