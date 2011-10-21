//=============================================================================
//
//  CLASS NConstraintGmmInterface
//
//=============================================================================


#ifndef COMISO_NCONSTRAINTGMMINTERFACE_HH
#define COMISO_NCONSTRAINTGMMINTERFACE_HH


//== INCLUDES =================================================================

#include <gmm/gmm.h>

#include <CoMISo/Config/CoMISoDefines.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

	      

/** \class NProblemGmmInterface NProblemGmmInterface.hh <ACG/.../NPRoblemGmmInterface.hh>

    Brief Description.
  
    A more elaborate description follows.
*/
class COMISODLLEXPORT NConstraintGmmInterface
{
public:
  
  // ToDo: appropriate Vector/MatrixType ???
  typedef gmm::wsvector<double>        SVectorNP;
  typedef gmm::row_matrix< SVectorNP > SMatrixNP;

  // different types of constraints
  enum ConstraintType {NC_EQUAL, NC_LESS_EQUAL, NC_GREATER_EQUAL};

  /// Default constructor
  NConstraintGmmInterface(const ConstraintType _type = NC_EQUAL) : type_(_type) {}
 
  /// Destructor
  ~NConstraintGmmInterface() {}

  virtual int             n_unknowns      (                                      ) = 0;
  virtual double          eval_constraint ( const double* _x                     ) = 0;
  virtual void            eval_gradient   ( const double* _x, SVectorNP& _g      ) = 0;
  virtual void            eval_hessian    ( const double* _x, SMatrixNP& _h      ) = 0;

  virtual ConstraintType  constraint_type (                                      ) { return type_; }

  virtual bool            is_satisfied    ( const double* _x, double _eps = 1e-6 )
  {
    switch( type_)
    {
      case NC_EQUAL        : return (fabs(eval_constraint(_x)) <=  _eps); break;
      case NC_LESS_EQUAL   : return (     eval_constraint(_x)  <=  _eps); break;
      case NC_GREATER_EQUAL: return (     eval_constraint(_x)  >= -_eps); break;
    }
    return false;
  }

  virtual double gradient_update_factor( const double* _x, double _eps = 1e-6 )
  {
    double val = eval_constraint(_x);
    bool   upper_bound_ok = ( val <=  _eps);
    bool   lower_bound_ok = ( val >= -_eps);

    if(upper_bound_ok)
    {
      if(lower_bound_ok || type_ == NC_LESS_EQUAL) return 0.0;
      else                                         return 1.0;
    }
    else
    {
      if(lower_bound_ok && type_ == NC_GREATER_EQUAL) return  0.0;
      else                                            return -1.0;
    }
  }


private:
  // constraint type
  ConstraintType type_;
};


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // ACG_NCONSTRAINTGMMINTERFACE_HH defined
//=============================================================================

