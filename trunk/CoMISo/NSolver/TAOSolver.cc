//=============================================================================
//
//  CLASS TAOSolver - IMPLEMENTATION
//
//=============================================================================

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_TAO_AVAILABLE


//== INCLUDES =================================================================

#include "TAOSolver.hh"

//#include <dlfcn.h>

//== NAMESPACES ===============================================================

namespace ACG {

//== IMPLEMENTATION ========================================================== 

// static member initialization
bool TAOSolver::initialized_ = false;



int
TAOSolver::
solve( NProblemGmmInterface* _base)
{
  // // initialize (only once)
  // initialize();

  // std::cerr << "tao 1\n";
  // //  MPI_Init(0,0); 
  // char *libm="libmpi.so"; 
  // dlopen(libm,RTLD_GLOBAL);

  if(!initialized_)
  {
    /* Initialize TAO,PETSc */
    // non command line arguments necessary ...
    std::cerr << "Initialize MPI/Petsc/TAO ";
    static char  help[] ="help\n";
    int argc = 0;
    char **argv;
    //    MPI_Init(&argc, &argv);
    PetscInitialize( &argc, &argv,(char *)0,help );
    TaoInitialize  ( &argc, &argv,(char *)0,help );
    
    initialized_ = true;
    std::cerr << " done!!!\n";
  }

  /* used to check for functions returning nonzeros */
  int             info;

  // check for single processor
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if (size >1) {
    PetscPrintf(PETSC_COMM_SELF,"TAOSolver is intended for single processor use!\n");
    SETERRQ(1,"Incorrect number of processors");
  }

  /* Create TAO solver and set desired solution method  */
  //  TaoMethod       method="tao_cg";  /* minimization method */
  TaoMethod       method="tao_ntr";  /* minimization method */
  //  TaoMethod       method="tao_nm";  /* minimization method */
  TAO_SOLVER      tao;               /* TAO_SOLVER solver context */
  TAO_APPLICATION testapp;        /* The PETSc application */

  info = TaoCreate(PETSC_COMM_SELF,method,&tao); CHKERRQ(info);
  info = TaoApplicationCreate(PETSC_COMM_SELF,&testapp); CHKERRQ(info);

  // initalize vector
  int n = _base->n_unknowns();
  Vec x;
  info = VecCreateSeq(PETSC_COMM_SELF, n, &x); CHKERRQ(info);
  PetscScalar* X;
  info = VecGetArray(x,&X); CHKERRQ(info);
  _base->initial_x(X);
  info = VecRestoreArray(x,&X); CHKERRQ(info);

  // initialize matrix
 /* Create a matrix data structure to store the Hessian.  This structure will be used by TAO */
  Mat H;
  // ToDo: get nonzero_pattern
  //  int nnz[1]; nnz[0] = 1;
  //  info = MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,0,nnz,&H);        /* PETSc routine */
  info = MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,0,0,&H);        /* PETSc routine */
  info = MatSetOption(H,MAT_SYMMETRIC,PETSC_TRUE); CHKERRQ(info); /* PETSc flag    */
  info = TaoAppSetHessianMat(testapp,H,H); CHKERRQ(info); /* A TAO routine */

  // initialize solution vector
  info = TaoAppSetInitialSolutionVec(testapp,x); CHKERRQ(info);

  /* Provide TAO routines for function evaluation */
  info = TaoAppSetObjectiveRoutine(testapp, objective, (void*) _base); CHKERRQ(info);
  info = TaoAppSetGradientRoutine (testapp, gradient , (void*) _base); CHKERRQ(info);
  info = TaoAppSetHessianRoutine  (testapp, hessian  , (void*) _base); CHKERRQ(info);

  /* SOLVE THE APPLICATION */
  info = TaoSolveApplication(testapp,tao); CHKERRQ(info);

  /* Get information on termination */
  TaoTerminateReason reason;
  info = TaoGetTerminationReason(tao,&reason); CHKERRQ(info);
  if (reason <= 0)
    std::cerr << "Warning: TAO-Solver did not converge!!!\n";
  else
    std::cerr << "TAO-Solver converged!!!\n";

  // To View TAO solver information use
  info = TaoView(tao); CHKERRQ(info); 

  // if successfull get and store result
  //  if( reason)
  {
    info = TaoAppGetSolutionVec(testapp, &x);
    info = VecGetArray(x,&X); CHKERRQ(info);
    _base->store_result( X);
    info = VecRestoreArray(x,&X); CHKERRQ(info);
  }
  //  VecView(x, PETSC_VIEWER_STDOUT_WORLD);

  // /* Free TAO data structures */
  info = TaoDestroy(tao); CHKERRQ(info);
  info = TaoAppDestroy(testapp); CHKERRQ(info);

  return reason;
}


//-----------------------------------------------------------------------------


int
TAOSolver::
objective( TAO_APPLICATION _app, Vec _x, double* _result, void* _base)
{
  NProblemGmmInterface* base = (NProblemGmmInterface*) _base;
  
  PetscScalar *x;

  /* Get pointers to vector data */
  int info = VecGetArray(_x,&x); CHKERRQ(info);

  // evaluate function
  (*_result) = base->eval_f(x);

  /* Restore vectors */
  info = VecRestoreArray(_x,&x); CHKERRQ(info);

  return 0;
}


//-----------------------------------------------------------------------------


int
TAOSolver::
gradient(TAO_APPLICATION _app, Vec _x, Vec _g, void* _base)
{
  NProblemGmmInterface* base = (NProblemGmmInterface*) _base;

  PetscScalar *x, *g;
  int info;

  /* Get pointers to vector data */
  info = VecGetArray(_x,&x); CHKERRQ(info);
  info = VecGetArray(_g,&g); CHKERRQ(info);

  // compute gradient
  base->eval_gradient( x, g);

  /* Restore vectors */
  info = VecRestoreArray(_x,&x); CHKERRQ(info);
  info = VecRestoreArray(_g,&g); CHKERRQ(info);

  return 0;
}


//-----------------------------------------------------------------------------


int
TAOSolver::
hessian(TAO_APPLICATION _app, Vec _x, Mat* _H, Mat* _H_pre, MatStructure* _H_struct, void* _base)
{
  NProblemGmmInterface* base = (NProblemGmmInterface*) _base;

  PetscScalar *x;

  /* Get pointers to vector data */
  int info = VecGetArray(_x,&x); CHKERRQ(info);

  /* Initialize matrix entries to zero */
  info = MatZeroEntries(*_H);  CHKERRQ(info);

  // iterate over non-zero elements
  NProblemGmmInterface::SMatrixNS H;
  base->eval_hessian( x, H);

  for (unsigned int i = 0; i < gmm::mat_nrows(H); ++i) 
  {
    typedef gmm::linalg_traits<NProblemGmmInterface::SMatrixNS>::const_sub_row_type
      CRow;
    CRow row = gmm::mat_const_row(H, i);

    gmm::linalg_traits<CRow>::const_iterator it  = gmm::vect_const_begin(row);
    gmm::linalg_traits<CRow>::const_iterator ite = gmm::vect_const_end(row);
    
    int m = 1;
    int n = 1;
    int idxm[1]; idxm[0] = i;
    int idxn[1]; 
    PetscScalar values[1]; 
    for(; it != ite; ++it)
    {
      idxn[0] = it.index();
      values[0] = *it;
      info = MatSetValues(*_H, m, idxm, n, idxn, values, INSERT_VALUES);
    }
  }

  /* Assemble the matrix */
  info = MatAssemblyBegin(*_H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);
  info = MatAssemblyEnd(*_H,MAT_FINAL_ASSEMBLY); CHKERRQ(info);

  *_H_struct = SAME_NONZERO_PATTERN;

  /* Restore vectors */
  info = VecRestoreArray(_x,&x); CHKERRQ(info);

  return 0;
}


//-----------------------------------------------------------------------------


void
TAOSolver::
initialize()
{
  if(!initialized_)
  {
    /* Initialize TAO,PETSc */
    // non command line arguments necessary ...
    std::cerr << "Initialize MPI/Petsc/TAO ";
    static char  help[] ="help\n";
    static int argc = 0;
    static char **argv;
    //    MPI_Init(&argc, &argv);
    PetscInitialize( &argc, &argv,(char *)0,help );
    TaoInitialize  ( &argc, &argv,(char *)0,help );
    
    initialized_ = true;
    std::cerr << " done!!!\n";
  }
}


//-----------------------------------------------------------------------------


void
TAOSolver::
cleanup()
{
  /* Finalize TAO */
  TaoFinalize();
  PetscFinalize();
}


//=============================================================================
} // namespace ACG
//=============================================================================
#endif // COMISO_TAO_AVAILABLE

