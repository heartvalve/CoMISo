//=============================================================================
//
//  CLASS TAOSolver
//
//=============================================================================


#ifndef COMISO_TAOSOLVER_HH
#define COMISO_TAOSOLVER_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_TAO_AVAILABLE

//== INCLUDES =================================================================

#include <mpi.h>
#include <tao.h>
#include <gmm/gmm.h>

#include "NProblemGmmInterface.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class TAOSolver TAOSolver.hh <ACG/.../TAOSolver.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class TAOSolver
{
public:
   
  /// Default constructor
  TAOSolver() {}
 
  /// Destructor
  ~TAOSolver() {}

  // solve problem
  static int solve( NProblemGmmInterface* _base);


private:

  static void initialize();

  // ToDo: cleanup has to be started automatically
  static void cleanup();

  // declar TAO function prototypes
  static int objective(TAO_APPLICATION,Vec,double*,void*);
  static int gradient (TAO_APPLICATION,Vec,Vec    ,void*);
  static int hessian  (TAO_APPLICATION,Vec,Mat*,Mat*,MatStructure*,void*);

private:
  // initialized?
  static bool initialized_;
};


//=============================================================================
} // namespace COMISO

//=============================================================================
#endif // COMISO_TAO_AVAILABLE
//=============================================================================
#endif // ACG_TAOSOLVER_HH defined
//=============================================================================

