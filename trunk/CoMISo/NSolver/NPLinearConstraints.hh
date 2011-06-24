//=============================================================================
//
//  CLASS NPLinearConstraints
//
//=============================================================================


#ifndef ACG_NPLINEARCONSTRAINTS_HH
#define ACG_NPLINEARCONSTRAINTS_HH


//== INCLUDES =================================================================

#include <gmm/gmm.h>
#include "NSolverGmmInterface.hh"
#include "LinearConstraintHandlerElimination.hh"
#include "LinearConstraintHandlerPenalty.hh"
//#include "LinearConstraintHandlerLagrange.hh"


//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================

	      

/** \class NPLinearConstraints{ Elimination, Penalty, Lagrange}
providing different techniques to handle linear constraints within
non-linear optimization. All three techniques transfer the constraint
problem into an unconstrained one.
*/


// ConstraintHandler working with Elimination approach
class NPLinearConstraintsElimination : public NSolverGmmInterface, public LinearConstraintHandlerElimination
{
public:

  // use c-arrays as vectors for gmm
  typedef gmm::array1D_reference<double*> VectorPT;

  /// Default constructor having a pointer to the main problem
  NPLinearConstraintsElimination( NSolverGmmInterface* _np) : base_(_np), cache_initialized_(false)
  {
    if( !base_) std::cerr << "Warning: initialized NPLinearConstraints with zero pointer...\n";
    this->initialize_identity(base_->n_unknowns());
  }
 
  /// Destructor
  ~NPLinearConstraintsElimination() {}

  // initialize constraints
  template<class MatrixT, class VectorT>
  void initialize_constraints( const MatrixT& _C, const VectorT& _c)
  {
    initialize(_C, _c);
    cache_initialized_ = false;
  }


  // NSolverInterface
  virtual int    n_unknowns   ()
  { return this->n_reduced(); }

  virtual void   initial_x    ( double* _x )
  {
    // transform initial x from parent NP
    x_.resize(this->n());
    base_->initial_x(P(x_));
    this->transform_x(P(x_), _x);
    cache_initialized_ = false;
    update_x(_x);
  }

  virtual double eval_f( const double* _x )
  {
    update_x(_x);

    return base_->eval_f(P(x_));
  }

  virtual void   eval_gradient( const double* _x, double*    _g)
  {
    update_x(_x);
    vtemp_.resize(this->n());
    base_->eval_gradient( P(x_), P(vtemp_));
    this->transform_gradient( P(vtemp_), _g);
  }

  virtual void   eval_hessian ( const double* _x, SMatrixNS& _H)
  {
    update_x(_x);
    SMatrixNS H;
    base_->eval_hessian(P(x_), H);
    this->transform_hessian(H,_H);
  }

  virtual void   store_result ( const double* _x )
  {
    update_x(_x);

    base_->store_result( P(x_));
  }

protected:

  void update_x(const double* _xr)
  {
    if(!cache_initialized_)
    {
      x_.resize(this->n());
      xr_.resize(this->n_reduced());

      if(!xr_.empty())
        xr_[0] = _xr[0] + 1.0;
      else
      {
        x_.resize(this->n());
        this->inv_transform_x(xr_, x_);
        cache_initialized_ = true;
      }
    }

    for( int i=0; i<this->n_reduced(); ++i)
      if( _xr[i] != xr_[i])
      {
        gmm::copy( VectorPT((double*)_xr, this->n_reduced()), this->xr_);
        x_.resize(this->n());
        this->inv_transform_x(xr_, x_);

        cache_initialized_ = true;

        //quit
        return;
      }
  }

  // convert vector into pointer
  double* P(const std::vector<double>& _v)
  {
    if(_v.size())
      return (double*) (&(_v[0]));
    else
      return 0;
  }

private:
  NSolverGmmInterface* base_;

  bool cache_initialized_;

  std::vector<double> x_;
  std::vector<double> xr_;
  std::vector<double> vtemp_;
};


//=============================================================================


// ConstraintHandler working with Penalty
class NPLinearConstraintsPenalty : public NSolverGmmInterface, public LinearConstraintHandlerPenalty
{
public:

  // use c-arrays as vectors for gmm
  typedef gmm::array1D_reference<double*> VectorPT;

  /// Default constructor having a pointer to the main problem
  NPLinearConstraintsPenalty( NSolverGmmInterface* _np) : base_(_np)
  { if( !base_) std::cerr << "Warning: initialized NPLinearConstraints with zero pointer...\n"; }

  /// Destructor
  ~NPLinearConstraintsPenalty() {}

  // initialize constraints
  template<class MatrixT, class VectorT>
  void initialize_constraints( const MatrixT& _C, const VectorT& _c)
  {
    initialize(_C, _c);
  }

  // NSolverInterface
  virtual int    n_unknowns   ()
  { return base_->n_unknowns(); }

  virtual void   initial_x    ( double* _x )
  {
    base_->initial_x(_x);
  }

  virtual double eval_f( const double* _x )
  {
    return this->add_penalty_f((double*)_x, base_->eval_f(_x));
  }

  virtual void   eval_gradient( const double* _x, double*    _g)
  {
    base_->eval_gradient( _x, _g);
    this->add_penalty_gradient((double*)_x, _g);
  }

  virtual void   eval_hessian ( const double* _x, SMatrixNS& _H)
  {
    base_->eval_hessian(_x, _H);
    this->add_penalty_hessian(_H);
  }

  virtual void   store_result ( const double* _x )
  {
    base_->store_result( _x);
  }

private:
  NSolverGmmInterface* base_;
};


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_NSOLVERGMMINTERFACE_HH defined
//=============================================================================

