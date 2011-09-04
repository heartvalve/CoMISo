//=============================================================================
//
//  CLASS IPOPTSolver - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_IPOPT_AVAILABLE
//=============================================================================


#include "IPOPTSolver.hh"

//== NAMESPACES ===============================================================

namespace COMISO {

//== IMPLEMENTATION ========================================================== 



int
IPOPTSolver::
solve(NProblemGmmInterface* _problem, std::vector<NConstraintGmmInterface*>& _constraints)
{
  //----------------------------------------------------------------------------
  // 1. Create an instance of IPOPT NLP
  //----------------------------------------------------------------------------
  Ipopt::SmartPtr<Ipopt::TNLP> np = new NProblemIPOPT(_problem, _constraints);

  //----------------------------------------------------------------------------
  // 2. solve problem
  //----------------------------------------------------------------------------
  // Create an instance of the IpoptApplication
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

  app->Options()->SetStringValue("linear_solver", "ma57");
  //  app->Options()->SetStringValue("derivative_test", "second-order");
  //  app->Options()->SetIntegerValue("print_level", 0);

  // Initialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Ipopt::Solve_Succeeded)
  {
    printf("\n\n*** Error IPOPT during initialization!\n");
  }

  //----------------------------------------------------------------------------
  // 3. solve problem
  //----------------------------------------------------------------------------
  status = app->OptimizeTNLP(np);

  //----------------------------------------------------------------------------
  // 4. output statistics
  //----------------------------------------------------------------------------
  if (status == Ipopt::Solve_Succeeded || status == Ipopt::Solved_To_Acceptable_Level)
  {
    // Retrieve some statistics about the solve
    Ipopt::Index iter_count = app->Statistics()->IterationCount();
    printf("\n\n*** IPOPT: The problem solved in %d iterations!\n", iter_count);

    Ipopt::Number final_obj = app->Statistics()->FinalObjective();
    printf("\n\n*** IPOPT: The final value of the objective function is %e.\n", final_obj);
  }

  return status;
}


//== IMPLEMENTATION PROBLEM INSTANCE==========================================================


bool NProblemIPOPT::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // number of variables
  n = problem_->n_unknowns();

  // number of constraints
  m = constraints_.size();

  // get nonzero structure
  std::vector<double> x(n);
  problem_->initial_x(&(x[0]));
  // ToDo: perturb x

  // nonzeros in the jacobian of C_ and the hessian of the lagrangian
  SVectorNP g;
  SMatrixNP H;
  problem_->eval_hessian(&(x[0]), H);
  nnz_jac_g = 0;
  nnz_h_lag = gmm::nnz(H);

  // clear old data
  jac_g_iRow_.clear();
  jac_g_jCol_.clear();
  h_lag_iRow_.clear();
  h_lag_jCol_.clear();

  // get non-zero structure of initial hessian
  // iterate over rows
  for( int i=0; i<n; ++i)
  {
    SVectorNP& ri = H.row(i);

    SVectorNP_citer v_it  = gmm::vect_const_begin(ri);
    SVectorNP_citer v_end = gmm::vect_const_end  (ri);

    for(; v_it != v_end; ++v_it)
    {
      // store lower triangular part only
      if( i >= (int)v_it.index())
      {
        h_lag_iRow_.push_back(i);
        h_lag_jCol_.push_back(v_it.index());
      }
    }
  }


  // get nonzero structure of constraints
  for( int i=0; i<m; ++i)
  {
    constraints_[i]->eval_gradient(&(x[0]),g);
    constraints_[i]->eval_hessian (&(x[0]),H);
    nnz_jac_g += gmm::nnz(g);
    nnz_h_lag += gmm::nnz(H);

    // build table
    SVectorNP_citer v_it  = gmm::vect_const_begin(g);
    SVectorNP_citer v_end = gmm::vect_const_end  (g);

    for(; v_it != v_end; ++v_it)
    {
      jac_g_iRow_.push_back(i);
      jac_g_jCol_.push_back(v_it.index());
    }

    for( int i=0; i<n; ++i)
    {
      SVectorNP& ri = H.row(i);

      v_it  = gmm::vect_const_begin(ri);
      v_end = gmm::vect_const_end  (ri);

      for(; v_it != v_end; ++v_it)
      {
        // store lower triangular part only
        if( i >= (int)v_it.index())
        {
          h_lag_iRow_.push_back(i);
          h_lag_jCol_.push_back(v_it.index());
        }
      }
    }
  }

  // store for error checking...
  nnz_jac_g_ = nnz_jac_g;
  nnz_h_lag_ = nnz_h_lag;

  // We use the standard fortran index style for row/col entries
  index_style = C_STYLE;

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  // first clear all variable bounds
  for( int i=0; i<n; ++i)
  {
    // x_l[i] = Ipopt::nlp_lower_bound_inf;
    // x_u[i] = Ipopt::nlp_upper_bound_inf;

    x_l[i] = -1.0e19;
    x_u[i] =  1.0e19;
  }

  // set bounds for constraints
  for( int i=0; i<m; ++i)
  {
    // enum ConstraintType {NC_EQUAL, NC_LESS_EQUAL, NC_GREATER_EQUAL};
    switch(constraints_[i]->constraint_type())
    {
      case NConstraintGmmInterface::NC_EQUAL         : g_u[i] = 0.0   ; g_l[i] =  0.0   ; break;
      case NConstraintGmmInterface::NC_LESS_EQUAL    : g_u[i] = 0.0   ; g_l[i] = -1.0e19; break;
      case NConstraintGmmInterface::NC_GREATER_EQUAL : g_u[i] = 1.0e19; g_l[i] =  0.0   ; break;
      default                                        :  g_u[i] = 1.0e19; g_l[i] = -1.0e19; break;
    }
  }

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  // get initial value of problem instance
  problem_->initial_x(x);

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // return the value of the objective function
  obj_value = problem_->eval_f(x);
  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  problem_->eval_gradient(x, grad_f);

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  // evaluate all constraint functions
  for( int i=0; i<m; ++i)
    g[i] = constraints_[i]->eval_constraint(x);

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  if (values == NULL)
  {
    // return the (cached) structure of the jacobian of the constraints
    gmm::copy(jac_g_iRow_, VectorPTi(iRow, jac_g_iRow_.size()));
    gmm::copy(jac_g_jCol_, VectorPTi(jCol, jac_g_jCol_.size()));
  }
  else
  {
    // return the values of the jacobian of the constraints

    // return the structure of the jacobian of the constraints
    // global index
    int gi = 0;
    SVectorNP g;

    for( int i=0; i<m; ++i)
    {
      constraints_[i]->eval_gradient(x, g);

      // iterate over non-zero values
      SVectorNP_citer it  = gmm::vect_const_begin(g);
      SVectorNP_citer ite = gmm::vect_const_end(g);

      for (; it != ite; ++it)
      {
        if(gi < nele_jac)
          values[gi] = *it;
        ++gi;
      }
    }

    if( gi != nele_jac)
      std::cerr << "Warning: number of non-zeros in Jacobian of C is incorrect: "
                << gi << " vs " << nele_jac << std::endl;
  }

  return true;
}


//-----------------------------------------------------------------------------


bool NProblemIPOPT::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
  if (values == NULL)
  {
    // return the (cached) structure of the hessian
    gmm::copy(h_lag_iRow_, VectorPTi(iRow, h_lag_iRow_.size()));
    gmm::copy(h_lag_jCol_, VectorPTi(jCol, h_lag_jCol_.size()));
  }
  else
  {
    // return values.

    // global index
    int gi = 0;

    // get hessian of problem
    SMatrixNP H;
    problem_->eval_hessian(x, H);

    for( int i=0; i<n; ++i)
    {
      SVectorNP& ri = H.row(i);

      SVectorNP_citer v_it  = gmm::vect_const_begin(ri);
      SVectorNP_citer v_end = gmm::vect_const_end  (ri);

      for(; v_it != v_end; ++v_it)
      {
        // store lower triangular part only
        if( i >= (int)v_it.index())
        {
          if( gi < nele_hess)
            values[gi] = obj_factor*(*v_it);
          ++gi;
        }
      }
    }

    // Hessians of Constraints
    for(unsigned int j=0; j<constraints_.size(); ++j)
    {
      constraints_[j]->eval_hessian(x, H);

      for( int i=0; i<n; ++i)
      {
        SVectorNP& ri = H.row(i);

        SVectorNP_citer v_it  = gmm::vect_const_begin(ri);
        SVectorNP_citer v_end = gmm::vect_const_end  (ri);

        for(; v_it != v_end; ++v_it)
        {
          // store lower triangular part only
          if( i >= (int)v_it.index())
          {
            if( gi < nele_hess)
              values[gi] = lambda[j]*(*v_it);
            ++gi;
          }
        }
      }
    }

    // error check
    if( gi != nele_hess)
      std::cerr << "Warning: number of non-zeros in Hessian of Lagrangian is incorrect: "
                << gi << " vs " << nele_hess << std::endl;
  }
  return true;
}


//-----------------------------------------------------------------------------


void NProblemIPOPT::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
                              const IpoptData* ip_data,
                              IpoptCalculatedQuantities* ip_cq)
{
  // problem knows what to do
  problem_->store_result(x);
}



//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_IPOPT_AVAILABLE
//=============================================================================