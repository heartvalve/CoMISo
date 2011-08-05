//=============================================================================
//
//  CLASS NPDERIVATIVECHECKER
//
//=============================================================================


#ifndef COMISO_NPDERIVATIVECHECKER_HH
#define COMISO_NPDERIVATIVECHECKER_HH


//== INCLUDES =================================================================

#include <iostream>
#include <iomanip>

#include <ACG/Utils/StopWatch.hh>
#include <gmm/gmm.h>
#include "NProblemGmmInterface.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class NPDerivativeChecker

    Brief Description.
  
    A more elaborate description follows.
*/
class NPDerivativeChecker
{
public:

  struct Config
  {
    Config() : x_min(-1.0), x_max(1.0), n_iters(1), dx(1e-5), eps(1e-3)
    {}

    double x_min;
    double x_max;
    int    n_iters;
    double dx;
    double eps;
  };
  
  /// Default constructor
  NPDerivativeChecker() {}

  /// Destructor
  ~NPDerivativeChecker() {}

  bool check_all(NProblemGmmInterface* _np, double _dx, double _eps)
  {
    conf_.dx = _dx;
    conf_.eps = _eps;
    return check_all(_np);
  }

  bool check_all(NProblemGmmInterface* _np)
  {
    bool d1_ok = check_d1(_np);
    bool d2_ok = check_d2(_np);

    return ( d1_ok && d2_ok);
  }

  bool check_d1(NProblemGmmInterface* _np)
  {
    int n_ok     = 0;
    int n_errors = 0;


    int n = _np->n_unknowns();
    std::vector<double> x(n), g(n);

    for(int i=0; i<conf_.n_iters; ++i)
    {
      // get random x
      get_random_x(x, conf_.x_min, conf_.x_max);
      // gradient
      _np->eval_gradient(P(x), P(g));

      for(int j=0; j<n; ++j)
      {
//        double fd = finite_difference(_np, x, j);
        x[j] += conf_.dx;
        double f1 = _np->eval_f(P(x));
        x[j] -= 2.0*conf_.dx;
        double f0 = _np->eval_f(P(x));
        x[j] += conf_.dx;
        double fd = (f1-f0)/(2.0*conf_.dx);

        if( fabs(fd-g[j]) > conf_.eps)
        {
          ++ n_errors;
          std::cerr << "Gradient error in component " << j << ": " << g[j]
                    << " should be " << fd << " (" << fabs(fd-g[j]) << ")" << std::endl;
        }
        else ++ n_ok;
      }
    }
    std::cerr << "############## Gradient Checker #############\n";
    std::cerr << "#ok   : " << n_ok << std::endl;
    std::cerr << "#error: " << n_errors << std::endl;

    return (n_errors == 0);
  }

  bool check_d2(NProblemGmmInterface* _np)
  {
    int n_ok     = 0;
    int n_errors = 0;

    int n = _np->n_unknowns();
    std::vector<double> x(n);
    NProblemGmmInterface::SMatrixNP H(n,n);

    for(int i=0; i<conf_.n_iters; ++i)
    {
      // get random x
      get_random_x(x, conf_.x_min, conf_.x_max);
      // gradient
      _np->eval_hessian(P(x), H);

      for(int j=0; j<n; ++j)
        for(int k=0; k<n; ++k)
        {
          x[j] += conf_.dx;
          x[k] += conf_.dx;
          double f0 = _np->eval_f(P(x));
          x[j] -= 2.0*conf_.dx;
          double f1 = _np->eval_f(P(x));
          x[j] += 2.0*conf_.dx;
          x[k] -= 2.0*conf_.dx;
          double f2 = _np->eval_f(P(x));
          x[j] -= 2.0*conf_.dx;
          double f3 = _np->eval_f(P(x));

          double fd = (f0-f1-f2+f3)/(4.0*conf_.dx*conf_.dx);


          if( fabs(fd-H(j,k)) > conf_.eps)
          {
            ++ n_errors;
            std::cerr << "Hessian error in component " << j << "," << k << ": " << H(j,k)
                      << " should be (following FiniteDifferences) " << fd << " (" << fabs(fd-H(j,k)) << ")" << std::endl;
          }
          else ++ n_ok;
        }
    }

    std::cerr << "############## Hessian Checker #############\n";
    std::cerr << "#ok   : " << n_ok << std::endl;
    std::cerr << "#error: " << n_errors << std::endl;

    return (n_errors == 0);
  }

  Config& config() { return conf_; }

protected:
  void get_random_x(std::vector<double>& _x, double _xmin, double _xmax)
  {
    // get random values in [-1,1]
    gmm::fill_random(_x);
    double range = _xmax - _xmin;
    for(unsigned int i=0; i<_x.size(); ++i)
      _x[i] = (((_x[i]+1.0)/2.0)*range + _xmin);
  }

  double* P(std::vector<double>& _v)
  {
    if( !_v.empty())
      return ((double*)&_v[0]);
    else
      return 0;
  }

private:
  Config conf_;
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // ACG_NPDERIVATIVECHECKER defined
//=============================================================================

