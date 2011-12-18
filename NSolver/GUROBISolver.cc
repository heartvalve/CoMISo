//=============================================================================
//
//  CLASS GUROBISolver - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_GUROBI_AVAILABLE
//=============================================================================


#include "GUROBISolver.hh"

//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION ========================================================== 


//// Constructor
//GUROBISolver::
//GUROBISolver()
//{
//
//}

//-----------------------------------------------------------------------------

// ********** SOLVE **************** //
bool
GUROBISolver::
solve(NProblemInterface*                  _problem,
      std::vector<NConstraintInterface*>& _constraints,
      std::vector<PairUiV>&               _discrete_constraints,
      const double                        _time_limit)
{
  try
  {
    //----------------------------------------------
    // 0. set up environment
    //----------------------------------------------

    GRBEnv   env   = GRBEnv();
    GRBModel model = GRBModel(env);

    model.getEnv().set(GRB_DoubleParam_TimeLimit, _time_limit);


    //----------------------------------------------
    // 1. allocate variables
    //----------------------------------------------

    // determine variable types: 0->real, 1->integer, 2->bool
    std::vector<char> vtypes(_problem->n_unknowns(),0);
    for(unsigned int i=0; i<_discrete_constraints.size(); ++i)
      switch(_discrete_constraints[i].second)
      {
        case Integer: vtypes[_discrete_constraints[i].first] = 1; break;
        case Binary : vtypes[_discrete_constraints[i].first] = 2; break;
        default     : break;
      }

    // GUROBI variables
    std::vector<GRBVar> vars;
    // first all
    for( int i=0; i<_problem->n_unknowns(); ++i)
      switch(vtypes[i])
      {
        case 0 : vars.push_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS) ); break;
        case 1 : vars.push_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_INTEGER   ) ); break;
        case 2 : vars.push_back( model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_BINARY    ) ); break;
      }


    // Integrate new variables
    model.update();

    //----------------------------------------------
    // 2. setup constraints
    //----------------------------------------------

    // get zero vector
    std::vector<double> x(_problem->n_unknowns(), 0.0);

    for(unsigned int i=0; i<_constraints.size();  ++i)
    {
      if(!_constraints[i]->is_linear())
        std::cerr << "Warning: GUROBISolver received a problem with non-linear constraints!!!" << std::endl;

      GRBLinExpr lin_expr;
      NConstraintInterface::SVectorNC gc;
      _constraints[i]->eval_gradient(P(x), gc);

      NConstraintInterface::SVectorNC::InnerIterator v_it(gc);
      for(; v_it; ++v_it)
//        lin_expr += v_it.value()*vars[v_it.index()];
        lin_expr = lin_expr + vars[v_it.index()]*v_it.value();

      double b = _constraints[i]->eval_constraint(P(x));

      switch(_constraints[i]->constraint_type())
      {
        case NConstraintInterface::NC_EQUAL         : model.addConstr(lin_expr + b == 0); break;
        case NConstraintInterface::NC_LESS_EQUAL    : model.addConstr(lin_expr + b <= 0); break;
        case NConstraintInterface::NC_GREATER_EQUAL : model.addConstr(lin_expr + b >= 0); break;
      }
    }
    model.update();

    //----------------------------------------------
    // 3. setup energy
    //----------------------------------------------

    if(!_problem->constant_hessian())
      std::cerr << "Warning: GUROBISolver received a problem with non-constant hessian!!!" << std::endl;

    GRBQuadExpr objective;

    // add quadratic part
    NProblemInterface::SMatrixNP H;
    _problem->eval_hessian(P(x), H);
    for( int i=0; i<H.outerSize(); ++i)
      for (NProblemInterface::SMatrixNP::InnerIterator it(H,i); it; ++it)
        objective += 0.5*it.value()*vars[it.row()]*vars[it.col()];


    // add linear part
    std::vector<double> g(_problem->n_unknowns());
    _problem->eval_gradient(P(x), P(g));
    for(unsigned int i=0; i<g.size(); ++i)
      objective += g[i]*vars[i];

    // add constant part
    objective += _problem->eval_f(P(x));

    model.set(GRB_IntAttr_ModelSense, 1);
    model.setObjective(objective);

    //----------------------------------------------
    // 4. solve problem
    //----------------------------------------------

    model.optimize();

    //----------------------------------------------
    // 5. store result
    //----------------------------------------------

    for(unsigned int i=0; i<vars.size(); ++i)
      x[i] = vars[i].get(GRB_DoubleAttr_X);

    _problem->store_result(P(x));

    std::cout << "GUROBI Objective: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    return true;
  }
  catch(GRBException e)
  {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    return false;
  }
  catch(...)
  {
    std::cout << "Exception during optimization" << std::endl;
    return false;
  }

  return false;
}


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_GUROBI_AVAILABLE
//=============================================================================
