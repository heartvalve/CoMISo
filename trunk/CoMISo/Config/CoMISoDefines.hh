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


#ifndef COMISODLLEXPORT
	#ifdef WIN32
		#ifdef COMISOMDLL
			#ifdef USECOMISO
				#define COMISODLLEXPORT __declspec(dllimport)
				#define COMISODLLEXPORTONLY 
			#else
				#define COMISODLLEXPORT __declspec(dllexport)
				#define COMISODLLEXPORTONLY __declspec(dllexport)
			#endif
		#else		
			#define COMISODLLEXPORT	
			#define COMISODLLEXPORTONLY
		#endif
	#else
		#define COMISODLLEXPORT
		#define COMISODLLEXPORTONLY
	#endif
#endif

#undef min
#undef max


