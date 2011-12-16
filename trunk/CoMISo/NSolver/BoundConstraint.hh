//=============================================================================
//
//  CLASS BoundConstraint
//
//=============================================================================


#ifndef COMISO_BOUNDCONSTRAINT_HH
#define COMISO_BOUNDCONSTRAINT_HH


//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include "NConstraintInterface.hh"


//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO { 

//== CLASS DEFINITION =========================================================

	      

/** \class BoundConstraint

    Brief Description.
  
    A more elaborate description follows.
*/
class COMISODLLEXPORT BoundConstraint : public NConstraintInterface
{
public:
  
//  inherited from NConstraintInterface
//  typedef Eigen::SparseVector<double> SVectorNC;
//  typedef SuperSparseMatrixT<double>  SMatrixNC;
//  // different types of constraints
//  enum ConstraintType {NC_EQUAL, NC_LESS_EQUAL, NC_GREATER_EQUAL};

  /// Default constructor
  BoundConstraint(const unsigned int  _var_idx = 0,   // index of variable for bound constraint
                  const double        _bound   = 0.0, // bound: x(_var_idx) #_type, <,=,># _bound
                  const unsigned int   _n      = 0,   // number of unknowns in problem
                  const ConstraintType _type   = NC_LESS_EQUAL) // type of bound upper, lower or both (equal)
  : NConstraintInterface(_type), idx_(_var_idx), bound_(_bound), n_(_n)
  {}
 
  /// Destructor
  ~BoundConstraint() {}

  virtual int    n_unknowns      (                                      ) { return n_;}
  virtual double eval_constraint ( const double* _x                     ) { return _x[idx_] - bound_; }
  virtual void   eval_gradient   ( const double* _x, SVectorNC& _g      ) { _g.resize(n_); _g.coeffRef(idx_) = 1.0; }
  virtual void   eval_hessian    ( const double* _x, SMatrixNC& _h      ) { _h.clear(); _h.resize(n_,n_); }


  // set/get values
  unsigned int& idx()   {return idx_;}
  double&       bound() {return bound_;}
  unsigned int& n()     {return n_;}

private:
  // variable idx
  unsigned int idx_;
  // variable bound
  double bound_;
  // number of unknowns
  unsigned int n_;
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // ACG_BOUNDCONSTRAINT_HH defined
//=============================================================================

