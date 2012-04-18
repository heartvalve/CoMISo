//=============================================================================
//
//  CLASS GUROBISolver - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include "CPLEXSolver.hh"
#if COMISO_CPLEX_AVAILABLE
//=============================================================================


#include <stdexcept>

//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION ========================================================== 



CPLEXSolver::
CPLEXSolver()
{
}

//-----------------------------------------------------------------------------

// ********** SOLVE **************** //
bool
CPLEXSolver::
solve(NProblemInterface*                  _problem,
      std::vector<NConstraintInterface*>& _constraints,
      std::vector<PairIndexVtype>&        _discrete_constraints,
      const double                        _time_limit)
{

  try
  {
    //----------------------------------------------
    // 0. set up environment
    //----------------------------------------------

    IloModel model(env_);
    //    model.getEnv().set(GRB_DoubleParam_TimeLimit, _time_limit);

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

    // CPLEX variables
    std::vector<IloNumVar> vars;
    // first all
    for( int i=0; i<_problem->n_unknowns(); ++i)
      switch(vtypes[i])
      {
        case 0 : vars.push_back( IloNumVar(env_,-IloInfinity, IloInfinity, IloNumVar::Float) ); break;
        case 1 : vars.push_back( IloNumVar(env_,  -IloIntMax,   IloIntMax, IloNumVar::Int)   ); break;
        case 2 : vars.push_back( IloNumVar(env_,           0,           1, IloNumVar::Bool)  ); break;
      }


    // Integrate new variables
//    model.update();

    //----------------------------------------------
    // 2. setup constraints
    //----------------------------------------------

    // get zero vector
    std::vector<double> x(_problem->n_unknowns(), 0.0);

    for(unsigned int i=0; i<_constraints.size();  ++i)
    {
      if(!_constraints[i]->is_linear())
        std::cerr << "Warning: GUROBISolver received a problem with non-linear constraints!!!" << std::endl;

      IloExpr lin_expr(env_);
      NConstraintInterface::SVectorNC gc;
      _constraints[i]->eval_gradient(P(x), gc);

      NConstraintInterface::SVectorNC::InnerIterator v_it(gc);
      for(; v_it; ++v_it)
//        lin_expr += IloNumExpr(vars[v_it.index()]);
        lin_expr += vars[v_it.index()]*v_it.value();

      double b = _constraints[i]->eval_constraint(P(x));


      switch(_constraints[i]->constraint_type())
      {
        case NConstraintInterface::NC_EQUAL         : model.add(lin_expr + b == 0); break;
        case NConstraintInterface::NC_LESS_EQUAL    : model.add(lin_expr + b <= 0); break;
        case NConstraintInterface::NC_GREATER_EQUAL : model.add(lin_expr + b >= 0); break;
      }

    }

//    model.update();

    //----------------------------------------------
    // 3. setup energy
    //----------------------------------------------

    if(!_problem->constant_hessian())
      std::cerr << "Warning: CPLEXSolver received a problem with non-constant hessian!!!" << std::endl;

//    GRBQuadExpr objective;
    IloExpr objective(env_);
//
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

    model.add(IloMinimize(env_,objective));

//    model.set(GRB_IntAttr_ModelSense, 1);
//    model.setObjective(objective);
//    model.update();


    //----------------------------------------------
    // 4. solve problem
    //----------------------------------------------

    IloCplex cplex(model);
    cplex.setParam(IloCplex::TiLim, _time_limit);
    cplex.solve();


//    if (solution_input_path_.empty())
//    {
//      if (!problem_env_output_path_.empty())
//      {
//        std::cout << "Writing problem's environment into file \"" << problem_env_output_path_ << "\"." << std::endl;
//        model.getEnv().writeParams(problem_env_output_path_);
//      }
//      if (!problem_output_path_.empty())
//      {
//        std::cout << "Writing problem into file \"" << problem_output_path_ << "\"." << std::endl;
//        GurobiHelper::outputModelToMpsGz(model, problem_output_path_);
//      }
//
//      model.optimize();
//    }
//    else
//    {
//        std::cout << "Reading solution from file \"" << solution_input_path_ << "\"." << std::endl;
//    }
//
    //----------------------------------------------
    // 5. store result
    //----------------------------------------------

    // store computed result
    for(unsigned int i=0; i<vars.size(); ++i)
      x[i] = cplex.getValue(vars[i]);

    _problem->store_result(P(x));

/*
    if (solution_input_path_.empty())
    {
      // store computed result
      for(unsigned int i=0; i<vars.size(); ++i)
        x[i] = vars[i].get(GRB_DoubleAttr_X);
    }
*/
//    else
//    {
//        std::cout << "Loading stored solution from \"" << solution_input_path_ << "\"." << std::endl;
//        // store loaded result
//        const size_t oldSize = x.size();
//        x.clear();
//        GurobiHelper::readSolutionVectorFromSOL(x, solution_input_path_);
//        if (oldSize != x.size()) {
//            std::cerr << "oldSize != x.size() <=> " << oldSize << " != " << x.size() << std::endl;
//            throw std::runtime_error("Loaded solution vector doesn't have expected dimension.");
//        }
//    }
//
//    _problem->store_result(P(x));
//
//    // ObjVal is only available if the optimize was called.
//    if (solution_input_path_.empty())
//        std::cout << "GUROBI Objective: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
    return true;
  }
  catch (IloException& e)
  {
     cerr << "Concert exception caught: " << e << endl;
     return false;
  }
  catch (...)
  {
     cerr << "Unknown exception caught" << endl;
     return false;
  }

  return false;
}


//-----------------------------------------------------------------------------


//void
//CPLEXSolver::
//set_problem_output_path( const std::string &_problem_output_path)
//{
//  problem_output_path_ = _problem_output_path;
//}
//
//
////-----------------------------------------------------------------------------
//
//
//void
//CPLEXSolver::
//set_problem_env_output_path( const std::string &_problem_env_output_path)
//{
//  problem_env_output_path_ = _problem_env_output_path;
//}
//
//
////-----------------------------------------------------------------------------
//
//
//void
//CPLEXSolver::
//set_solution_input_path(const std::string &_solution_input_path)
//{
//  solution_input_path_ = _solution_input_path;
//}


//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_CPLEX_AVAILABLE
//=============================================================================
