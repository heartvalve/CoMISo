//=============================================================================
//
//  CLASS GUROBISolver
//
//=============================================================================


#ifndef COMISO_GUROBISOLVER_HH
#define COMISO_GUROBISOLVER_HH


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_GUROBI_AVAILABLE

//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include <vector>
#include "NProblemInterface.hh"
#include "NConstraintInterface.hh"
#include "VariableType.hh"

#include <gurobi_c++.h>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class NewtonSolver GUROBISolver.hh

    Brief Description.
  
    A more elaborate description follows.
*/
class COMISODLLEXPORT GUROBISolver
{
public:
   
  typedef std::pair<unsigned int, VariableType> PairUiV;

  /// Default constructor -> set up IpOptApplication
  GUROBISolver() {}
 
  /// Destructor
  ~GUROBISolver() {}

  // ********** SOLVE **************** //
  bool solve(NProblemInterface*                  _problem,                // problem instance
             std::vector<NConstraintInterface*>& _constraints,            // linear constraints
             std::vector<PairUiV>&               _discrete_constraints,   // discrete constraints
             const double                        _time_limit = 60     ); // time limit in seconds

protected:
  double* P(std::vector<double>& _v)
  {
    if( !_v.empty())
      return ((double*)&_v[0]);
    else
      return 0;
  }

private:

};



//=============================================================================
} // namespace COMISO

//=============================================================================
#endif // COMISO_GUROBI_AVAILABLE
//=============================================================================
#endif // ACG_GUROBISOLVER_HH defined
//=============================================================================

