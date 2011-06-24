//=============================================================================
//
//  CLASS BaseTaoGmmInterface
//
//=============================================================================


#ifndef ACG_NSOLVERGMMINTERFACE_HH
#define ACG_NSOLVERGMMINTERFACE_HH


//== INCLUDES =================================================================

#include <gmm/gmm.h>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================

	      

/** \class NSolverGmmInterface NSolverGmmInterface.hh <ACG/.../NSolverGmmInterface.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class NSolverGmmInterface
{
public:
  
  // ToDo: appropriate MatrixType ???
  typedef gmm::row_matrix< gmm::wsvector<double> > SMatrixNS;

  /// Default constructor
  NSolverGmmInterface() {}
 
  /// Destructor
  ~NSolverGmmInterface() {}

  virtual int    n_unknowns   (                                ) = 0;
  virtual void   initial_x    (       double* _x               ) = 0;
  virtual double eval_f       ( const double* _x               ) = 0;
  virtual void   eval_gradient( const double* _x, double*    _g) = 0;
  virtual void   eval_hessian ( const double* _x, SMatrixNS& _H) = 0;
  virtual void   store_result ( const double* _x               ) = 0;
};


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_NSOLVERGMMINTERFACE_HH defined
//=============================================================================

