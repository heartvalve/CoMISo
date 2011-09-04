//=============================================================================
//
//  CLASS NConstraintGmmInterface
//
//=============================================================================


#ifndef COMISO_LINEARCONSTRAINT_HH
#define COMISO_LINEARCONSTRAINT_HH


//== INCLUDES =================================================================

#include <gmm/gmm.h>
#include "NConstraintGmmInterface.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class NProblemGmmInterface NProblemGmmInterface.hh <ACG/.../NPRoblemGmmInterface.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class LinearConstraint : public NConstraintGmmInterface
{
public:
  
  // ToDo: appropriate Vector/MatrixType ???
  typedef gmm::wsvector<double>        SVectorNP;
  typedef gmm::row_matrix< SVectorNP > SMatrixNP;

  // use c-arrays as vectors for gmm
  typedef gmm::array1D_reference<double*> VectorPT;

  // different types of constraints
//  enum ConstraintType {NC_EQUAL, NC_LESS_EQUAL, NC_GREATER_EQUAL};

  /// Default constructor
  LinearConstraint(const ConstraintType _type = NC_EQUAL) : NConstraintGmmInterface(_type)
  {}

  // linear equation of the form -> coeffs_^T * (x,1) =_type= 0
  LinearConstraint(const SVectorNP& _coeffs_plus_const, const ConstraintType _type = NC_EQUAL) : NConstraintGmmInterface(_type)
  {
    int n = gmm::vect_size(_coeffs_plus_const)-1;

    gmm::resize(coeffs_, n+1);
    gmm::copy(_coeffs_plus_const, coeffs_);

    if( n >= 0)
      b_ = coeffs_[n];
    else
      b_ = 0.0;

    gmm::resize(coeffs_, n);
  }

  // linear equation of the form -> coeffs_^T * (x,1) =_type= 0
  LinearConstraint(const SVectorNP& _coeffs, const double _b, const ConstraintType _type = NC_EQUAL) : NConstraintGmmInterface(_type)
  {
    int n = gmm::vect_size(_coeffs);
    gmm::resize(coeffs_, n);
    gmm::copy(_coeffs, coeffs_);
    b_ = _b;
  }

  /// Destructor
  ~LinearConstraint() {}

  virtual int n_unknowns()
  {
    return gmm::vect_size(coeffs_);
  }

  virtual double eval_constraint ( const double* _x )
  {
     return (gmm::vect_sp(coeffs_, VectorPT((double*)_x, gmm::vect_size(coeffs_))) + b_);
  }

  virtual void eval_gradient( const double* _x, SVectorNP& _g      )
  {
    gmm::resize(_g, gmm::vect_size(coeffs_));
    gmm::copy  (coeffs_, _g);
  }

  virtual void eval_hessian    ( const double* _x, SMatrixNP& _h      )
  {
    gmm::resize(_h, gmm::vect_size(coeffs_), gmm::vect_size(coeffs_));
    gmm::clear(_h);
  }

  // inherited from base
//  virtual ConstraintType  constraint_type (                                      ) { return type_; }

private:

  // linear equation of the form -> coeffs_^T * x + b_
  SVectorNP coeffs_;
  double    b_;
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // ACG_LINEARCONSTRAINT_HH defined
//=============================================================================

