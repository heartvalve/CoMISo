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


//=============================================================================
//
//  CLASS CholmodSolver
//
//=============================================================================

#ifndef ACG_CHOLMOD_SOLVER_HH
#define ACG_CHOLMOD_SOLVER_HH


//== INCLUDES =================================================================


#include <CoMISo/Config/CoMISoDefines.hh>
#include "GMM_Tools.hh"

#include <iostream>
#include <vector>

typedef struct cholmod_common_struct cholmod_common;
typedef struct cholmod_factor_struct cholmod_factor;


//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================
class COMISODLLEXPORT CholmodSolver
{
public:

    // _size is maximal size this instance can handle (smaller problems are possible!!!)
    CholmodSolver();
    ~CholmodSolver();
    
    bool calc_system( const std::vector<int>&    _colptr, 
		      const std::vector<int>&    _rowind, 
		      const std::vector<double>& _values );


    template< class GMM_MatrixT>
    bool calc_system_gmm( const GMM_MatrixT& _mat);


    bool update_system( const std::vector<int>&    _colptr, 
 			const std::vector<int>&    _rowind, 
 			const std::vector<double>& _values );


    template< class GMM_MatrixT>
    bool update_system_gmm( const GMM_MatrixT& _mat);


    bool solve ( double *             _x0, double *             _b);

    bool solve ( std::vector<double>& _x0, std::vector<double>& _b)
    {return solve( &(_x0[0]), &(_b[0]));}

    
private:

    cholmod_common * mp_cholmodCommon;

    cholmod_factor * mp_L;

    std::vector<double> values_;
    std::vector<int>    colptr_;
    std::vector<int>    rowind_;

};

//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_CHOLMOD_SOLVER_TEMPLATES_C)
#define ACG_CHOLMOD_SOLVER_TEMPLATES
#include "CholmodSolverT.cc"
#endif
//=============================================================================
#endif // ACG_CHOLMOD_SOLVER_HH defined
//=============================================================================
