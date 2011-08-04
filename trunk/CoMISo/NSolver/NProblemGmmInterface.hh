//=============================================================================
//
//  CLASS NProblemGmmInterface
//
//=============================================================================


#ifndef ACG_NPROBLEMGMMINTERFACE_HH
#define ACG_NPROBLEMGMMINTERFACE_HH


//== INCLUDES =================================================================

#include <gmm/gmm.h>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================

	      

/** \class NProblemGmmInterface NProblemGmmInterface.hh <ACG/.../NPRoblemGmmInterface.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class NProblemGmmInterface
{
public:
  
  // ToDo: appropriate MatrixType ???
  typedef gmm::row_matrix< gmm::wsvector<double> > SMatrixNP;

  /// Default constructor
  NProblemGmmInterface() {}
 
  /// Destructor
  ~NProblemGmmInterface() {}

  virtual int    n_unknowns   (                                ) = 0;
  virtual void   initial_x    (       double* _x               ) = 0;
  virtual double eval_f       ( const double* _x               ) = 0;
  virtual void   eval_gradient( const double* _x, double*    _g) = 0;
  virtual void   eval_hessian ( const double* _x, SMatrixNP& _H) = 0;
  virtual void   store_result ( const double* _x               ) = 0;
};


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_NROBLEMGMMINTERFACE_HH defined
//=============================================================================

