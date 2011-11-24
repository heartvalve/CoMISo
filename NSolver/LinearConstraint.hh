//=============================================================================
//
//  CLASS NConstraintGmmInterface
//
//=============================================================================


#ifndef COMISO_LINEARCONSTRAINT_HH
#define COMISO_LINEARCONSTRAINT_HH


//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include "NConstraintInterface.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class NProblemGmmInterface NProblemGmmInterface.hh <ACG/.../NPRoblemGmmInterface.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class COMISODLLEXPORT LinearConstraint : public NConstraintInterface
{
public:

  // use c-arrays as vectors for gmm
  typedef gmm::array1D_reference<double*> VectorPT;
  typedef NConstraintInterface::SVectorNC SVectorNC;

  // different types of constraints
//  enum ConstraintType {NC_EQUAL, NC_LESS_EQUAL, NC_GREATER_EQUAL};

  /// Default constructor
  LinearConstraint(const ConstraintType _type = NC_EQUAL);

  // linear equation of the form -> coeffs_^T *x  + b_=_type= 0
  LinearConstraint(const SVectorNC& _coeffs, const double _b, const ConstraintType _type = NC_EQUAL);

  /// Destructor
  ~LinearConstraint();

  virtual int n_unknowns();

  const SVectorNC& coeffs() const;
        SVectorNC& coeffs();

  const double&    b() const;
        double&    b();

  virtual double eval_constraint ( const double* _x );
  
  virtual void eval_gradient( const double* _x, SVectorNC& _g      );

  virtual void eval_hessian    ( const double* _x, SMatrixNC& _h      );

  // inherited from base
//  virtual ConstraintType  constraint_type (                                      ) { return type_; }

private:

  // linear equation of the form -> coeffs_^T * x + b_
  SVectorNC coeffs_;
  double    b_;
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // ACG_LINEARCONSTRAINT_HH defined
//=============================================================================

