//=============================================================================
//
//  CLASS NewtonSolver
//
//=============================================================================


#ifndef ACG_NEWTONSOLVER_HH
#define ACG_NEWTONSOLVER_HH


//== INCLUDES =================================================================

#include <gmm/gmm.h>
#include "NSolverGmmInterface.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================

	      

/** \class NewtonSolver NewtonSolver.hh <ACG/.../NewtonSolver.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class NewtonSolver
{
public:
   
  /// Default constructor
  NewtonSolver() : max_iter_(20), convergence_eps_(1e-6), constant_hessian_structure_(false) {}
 
  /// Destructor
  ~NewtonSolver() {}

  // solve
  int solve(NSolverGmmInterface* _problem);

  // solve specifying parameters
  int solve(NSolverGmmInterface* _problem, int _max_iter, double _eps)
  {
    max_iter_ = _max_iter;
    convergence_eps_ = _eps;
    return solve(_problem);
  }

  bool& constant_hessian_structure() { return constant_hessian_structure_; }

protected:
  double* P(std::vector<double>& _v)
  {
    if( !_v.empty())
      return ((double*)&_v[0]);
    else
      return 0;
  }

private:
  int    max_iter_;
  double convergence_eps_;
  bool   constant_hessian_structure_;
};


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_NEWTONSOLVER_HH defined
//=============================================================================

