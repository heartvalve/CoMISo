/*===========================================================================*\
 *                                                                           *
 *                              CoMISo                                       *
 *      Copyright (C) 2008-2009 by Computer Graphics Group, RWTH Aachen      *
 *                      www.graphics.rwth-aachen.de                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is a part of CoMISo.                                          *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 1    $                                                       *
 *   $Author: zimmer $                                                       *
 *   $Date: 2009-08-05 16:35:37 +0200 (Wed, 05 Aug 2009) $                   *
 *                                                                           *
\*===========================================================================*/


#define ACG_CHOLMOD_SOLVER_TEMPLATES_C

#include "CholmodSolver.hh"


namespace ACG {


template< class GMM_MatrixT>
bool CholmodSolver::calc_system_gmm( const GMM_MatrixT& _mat)
{
//   std::vector<int>    colptr;
//   std::vector<int>    rowind;
//   std::vector<double> values;
    
    gmm::get_ccs_symmetric_data( _mat,
				 'l',
				 values_, 
				 rowind_, 
				 colptr_ );

    return calc_system( colptr_, rowind_, values_);
}
  

  //-----------------------------------------------------------------------------


template< class GMM_MatrixT>
bool CholmodSolver::update_system_gmm( const GMM_MatrixT& _mat)
{
//   std::vector<int>    colptr;
//   std::vector<int>    rowind;
//   std::vector<double> values;
    
  gmm::get_ccs_symmetric_data( _mat,
			       'l',
			       values_, 
			       rowind_, 
			       colptr_ );

    return update_system( colptr_, rowind_, values_);
}


}
