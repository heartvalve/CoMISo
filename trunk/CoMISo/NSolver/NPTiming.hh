//=============================================================================
//
//  CLASS NPTiming
//
//=============================================================================


#ifndef ACG_NPTIMING_HH
#define ACG_NPTIMING_HH


//== INCLUDES =================================================================

#include <iostream>
#include <iomanip>

#include <ACG/Utils/StopWatch.hh>
#include <gmm/gmm.h>
#include "NProblemGmmInterface.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace ACG {

//== CLASS DEFINITION =========================================================

	      

/** \class NProblemGmmInterface NProblemGmmInterface.hh <ACG/.../NProblemGmmInterface.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class NPTiming : public NProblemGmmInterface
{
public:
  
  /// Default constructor
  NPTiming(NProblemGmmInterface* _base) : base_(_base) {start_timing();}
 
  /// Destructor
  ~NPTiming() {}

  virtual int    n_unknowns   ()
  {
    return base_->n_unknowns();
  }

  virtual void   initial_x( double* _x )
  {
    base_->initial_x(_x);
  }

  virtual double eval_f( const double* _x )
  {
    ++n_eval_f_;
    sw_.start();
    double f = base_->eval_f(_x);
    timing_eval_f_ += sw_.stop();
    return f;
  }

  virtual void   eval_gradient( const double* _x, double*    _g)
  {
    ++n_eval_gradient_;
    sw_.start();
    base_->eval_gradient(_x, _g);
    timing_eval_gradient_ += sw_.stop();
  }

  virtual void   eval_hessian ( const double* _x, SMatrixNS& _H)
  {
    ++n_eval_hessian_;
    sw_.start();
    base_->eval_hessian(_x, _H);
    timing_eval_hessian_ += sw_.stop();
  }

  virtual void   store_result ( const double* _x )
  {
    base_->store_result(_x);
    print_statistics();
  }

  void start_timing()
  {
    swg_.start();

    timing_eval_f_ = 0.0;
    timing_eval_gradient_ = 0.0;
    timing_eval_hessian_ = 0.0;

    n_eval_f_ = 0;
    n_eval_gradient_ = 0;
    n_eval_hessian_ = 0;
  }

protected:

  void print_statistics()
  {
    double time_total = swg_.stop();

    double time_np = timing_eval_f_ + timing_eval_gradient_ + timing_eval_hessian_;



    std::cerr << "######## NP-Timings ########" << std::endl;
    std::cerr << "total time    : " << time_total/1000.0 << "s\n";
    std::cerr << "total time NP : " << time_np/1000.0 << "s  (" << time_np/time_total*100.0 << " %)\n";

    std::cerr << std::fixed << std::setprecision(5)
              << "eval_f time   : " << timing_eval_f_/1000.0
              << "s  ( #evals: " << n_eval_f_ << " -> avg "
              << timing_eval_f_/(1000.0*double(n_eval_f_)) << "s )\n"
              << "eval_grad time: " << timing_eval_gradient_/1000.0
              << "s  ( #evals: " << n_eval_gradient_ << " -> avg "
              << timing_eval_gradient_/(1000.0*double(n_eval_gradient_)) << "s )\n"
              << "eval_hess time: " << timing_eval_hessian_/1000.0
              << "s  ( #evals: " << n_eval_hessian_ << " -> avg "
              << timing_eval_hessian_/(1000.0*double(n_eval_hessian_)) << "s )\n";
  }

private:
  NProblemGmmInterface* base_;
  StopWatch swg_;
  StopWatch sw_;

  // timings
  double timing_eval_f_;
  double timing_eval_gradient_;
  double timing_eval_hessian_;

  // number of function executions
  int n_eval_f_;
  int n_eval_gradient_;
  int n_eval_hessian_;
};


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // ACG_NPTIMING_HH defined
//=============================================================================

