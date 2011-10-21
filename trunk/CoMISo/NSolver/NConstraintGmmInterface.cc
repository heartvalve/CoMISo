//=============================================================================
//
//  CLASS NConstraintGmmInterface
//
//=============================================================================


//== INCLUDES =================================================================

#include "NConstraintGmmInterface.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION =========================================================

	      
NConstraintGmmInterface::NConstraintGmmInterface(const ConstraintType _type ) :
type_(_type) {}

NConstraintGmmInterface::~NConstraintGmmInterface() {}

NConstraintGmmInterface::ConstraintType  NConstraintGmmInterface::constraint_type ( ) { return type_; }

bool NConstraintGmmInterface::is_satisfied    ( const double* _x, double _eps  )
{
  switch( type_)
  {
    case NC_EQUAL        : return (fabs(eval_constraint(_x)) <=  _eps); break;
    case NC_LESS_EQUAL   : return (     eval_constraint(_x)  <=  _eps); break;
    case NC_GREATER_EQUAL: return (     eval_constraint(_x)  >= -_eps); break;
  }
  return false;
}

double NConstraintGmmInterface::gradient_update_factor( const double* _x, double _eps )
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



//=============================================================================
} // namespace COMISO
//=============================================================================


