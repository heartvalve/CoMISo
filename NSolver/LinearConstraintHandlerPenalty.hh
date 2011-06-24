//=============================================================================
//
//  CLASS LinearConstraintHandlerElimination
//
//=============================================================================


#ifndef ACG_LINEARCONSTRAINTHANDLERPENALTY_HH
#define ACG_LINEARCONSTRAINTHANDLERPENALTY_HH


//== INCLUDES =================================================================

#include <iostream>
#include <gmm/gmm.h>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================



	      
/** \class LinearConstraintHandler LinearConstraintHandler.hh <ACG/.../LinearConstraintHandler.hh>

    Brief Description.
  
    A more elaborate description follows.
*/


class LinearConstraintHandlerPenalty
{
public:
   
  typedef gmm::col_matrix< gmm::wsvector< double> > CMatrix;
  typedef gmm::row_matrix< gmm::wsvector< double> > RMatrix;

  // use c-arrays as vectors for gmm
  typedef gmm::array1D_reference<double*> VectorPT;

  /// Constructor
  LinearConstraintHandlerPenalty() : penalty_(10000) {}

  // initialize Constructor
  template<class MatrixT, class VectorT>
  LinearConstraintHandlerPenalty( const MatrixT& _C, const VectorT& _c) : penalty_(10000)
  {initialize(_C, _c); }

 
  /// Destructor
  ~LinearConstraintHandlerPenalty() {}

  // penalty weight
  double& penalty() { return penalty_; }

  // number of variables
  int n() {return n_;}

  // number of linearly independent constraints (n-n_reduced)
  int n_constraints() { return m_;}

  // initialize new constraints
  template<class MatrixT, class VectorT>
  void initialize( const MatrixT& _C, const VectorT& _c);

  // initialize new constraints rhs only
  void initialize( const std::vector<double>& _c);
  void initialize( double* _c);

  // transform energy
  double add_penalty_f( double* _x, const double _f);

  // transform gradient
  void add_penalty_gradient( const std::vector<double>& _x, std::vector<double>& _g);
  void add_penalty_gradient( double* _x, double* _g);

  // transform hessian
  void add_penalty_hessian( RMatrix& _H);

private:
  // penalty weight
  double penalty_;
  
  // Linear Constraints C_*x_ = b_
  RMatrix             C_;
  std::vector<double> b_;

  // precomputed penalty terms
  RMatrix penalty_H_;
  std::vector<double> penalty_grad_b_;

  // temp vector
  std::vector<double> temp_;

  // number of variables
  int n_;
  // number of constraints (linear independent)
  int m_;
};

//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_LINEARCONSTRAINTHANDLERPENALTY_C)
#define ACG_LINEARCONSTRAINTHANDLERPENALTY_TEMPLATES
#include "LinearConstraintHandlerPenaltyT.cc"
#endif
//=============================================================================
#endif // ACG_LINEARCONSTRAINTHANDLERPENALTY_HH defined
//=============================================================================

