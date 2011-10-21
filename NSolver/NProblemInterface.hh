//=============================================================================
//
//  CLASS NProblemGmmInterface
//
//=============================================================================


#ifndef COMISO_NPROBLEMINTERFACE_HH
#define COMISO_NPROBLEMINTERFACE_HH


//== INCLUDES =================================================================

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

#include <CoMISo/Config/CoMISoDefines.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class NProblemGmmInterface NProblemGmmInterface.hh <ACG/.../NPRoblemGmmInterface.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class COMISODLLEXPORT NProblemInterface
{
public:
  
  // Sparse Matrix Type
  typedef Eigen::DynamicSparseMatrix<double,ColMajor> SMatrixNP;

  /// Default constructor
  NProblemInterface() {}
 
  /// Destructor
  ~NProblemInterface() {}

  virtual int    n_unknowns   (                                ) = 0;
  virtual void   initial_x    (       double* _x               ) = 0;
  virtual double eval_f       ( const double* _x               ) = 0;
  virtual void   eval_gradient( const double* _x, double*    _g) = 0;
  virtual void   eval_hessian ( const double* _x, SMatrixNP& _H) = 0;
  virtual void   store_result ( const double* _x               ) = 0;
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_NPROBLEMGMMINTERFACE_HH defined
//=============================================================================

