//=============================================================================
//
//  CLASS LinearProblem
//
//=============================================================================


#ifndef COMISO_LINEARPROBLEM_HH
#define COMISO_LINEARPROBLEM_HH


//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_EIGEN3_AVAILABLE

//== INCLUDES =================================================================

#include <stdio.h>
#include <iostream>
#include <vector>

#include <CoMISo/Config/CoMISoDefines.hh>
#include <CoMISo/NSolver/NProblemInterface.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class LinearProblem LinearProblem.hh <CoMISo/NSolver/LinearProblem.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class COMISODLLEXPORT LinearProblem : public NProblemInterface
{
public:
  
  /// Default constructor
  LinearProblem (unsigned int _dimension = 0)
  {
    // resize and zero elements
    coeffs_.resize(_dimension,0.0);
    x_.resize(_dimension, 0.0);
  }
 
  /// Destructor
  virtual ~LinearProblem() {}

  // problem definition
  virtual int    n_unknowns()
  { return coeffs_.size(); }

  virtual void   initial_x(double* _x)
  {
    if(!x_.empty())
      memcpy(_x, &(x_[0]), x_.size()*sizeof(double));
  }

  virtual double eval_f( const double* _x)
  {
    double d(0.0);
    for(unsigned int i=0; i<coeffs_.size(); ++i)
      d += coeffs_[i]*_x[i];
    return d;
  }

  virtual void   eval_gradient( const double* _x, double* _g)
  {
    if(!coeffs_.empty())
      memcpy(_g, &(coeffs_[0]), coeffs_.size()*sizeof(double));
  }

  virtual void   eval_hessian ( const double* _x, SMatrixNP& _H)
  {
    // resize and set to zero
    _H.resize(n_unknowns(), n_unknowns());
    _H.setZero();
  }

  virtual void   store_result ( const double* _x )
  {
    if(!x_.empty())
      memcpy(&(x_[0]), _x, x_.size()*sizeof(double));
  }

  // advanced properties
  virtual bool   constant_gradient() const { return true; }
  virtual bool   constant_hessian()  const { return true; }

  // give access to coefficients to define the problem
  std::vector<double>& coeffs() { return coeffs_;}

  // give access to result
  std::vector<double>& x() { return x_;}

private:

  // coefficients of linear problem
  std::vector<double> coeffs_;
  std::vector<double> x_;
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_EIGEN3_AVAILABLE
//=============================================================================
#endif // COMISO_NPROBLEMGMMINTERFACE_HH defined
//=============================================================================

