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



#ifndef VSTOOLS_HH
#define VSTOOLS_HH


//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

//== DEFINITION =========================================================

/** These functions are required for Visual Studio to work around missing 
    functions. Basic equivalent functions for double exist in the float 
    header but are named different. So this wrapper makes them standard compatible.
    */
#ifdef WIN32
 #include <float.h>

 namespace std {

   inline int isnan(double x)
   {
     return _isnan(x);
   } 

   inline int isinf(double x)
   {
     return !_finite(x);
   } 

  }

 inline double nearbyint(double x) {
   return double(int( x + 0.5));
 }

 inline double round ( double _value ) {
   return floor( _value + 0.5 );
 }


#endif


//=============================================================================
#endif // VSTOOLS_HH defined
//=============================================================================

