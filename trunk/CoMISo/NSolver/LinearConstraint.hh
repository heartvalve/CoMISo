//=============================================================================
//
//  CLASS NConstraintGmmInterface
//
//=============================================================================


#ifndef COMISO_LINEARCONSTRAINT_HH
#define COMISO_LINEARCONSTRAINT_HH


//== INCLUDES =================================================================

#include "NConstraintInterface.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class NProblemGmmInterface NProblemGmmInterface.hh <ACG/.../NPRoblemGmmInterface.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class LinearConstraint : public NConstraintInterface
{
public:

  // use c-arrays as vectors for gmm
  typedef gmm::array1D_reference<double*> VectorPT;

  // different types of constraints
//  enum ConstraintType {NC_EQUAL, NC_LESS_EQUAL, NC_GREATER_EQUAL};

  /// Default constructor
  LinearConstraint(const ConstraintType _type = NC_EQUAL) : NConstraintInterface(_type)
  {}

  // linear equation of the form -> coeffs_^T * (x,1) =_type= 0
  LinearConstraint(const SVectorNC& _coeffs, const double _b, const ConstraintType _type = NC_EQUAL) : NConstraintInterface(_type)
  {
    coeffs_ = _coeffs;
    b_ = _b;
  }

  /// Destructor
  ~LinearConstraint() {}

  virtual int n_unknowns()
  {
    return coeffs_.innerSize();
  }

  virtual double eval_constraint ( const double* _x )
  {
    double v = b_;

    SVectorNC::InnerIterator c_it(coeffs_);
    for(; c_it; ++c_it)
      v += c_it.value()*_x[c_it.index()];

    return v;
  }

  virtual void eval_gradient( const double* _x, SVectorNC& _g      )
  {
    _g = coeffs_;
  }

  virtual void eval_hessian    ( const double* _x, SMatrixNC& _h      )
  {
    _h.clear();
    _h.resize(coeffs_.innerSize(), coeffs_.innerSize());
  }

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

