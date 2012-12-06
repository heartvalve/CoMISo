/*
 * NProblemInterfaceAD.hpp
 *
 *  Created on: Nov 30, 2012
 *      Author: kremer
 */

#ifndef NPROBLEMINTERFACEAD_HPP_
#define NPROBLEMINTERFACEAD_HPP_

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_ADOLC_AVAILABLE
#if COMISO_Eigen3_AVAILABLE

//== INCLUDES =================================================================

#include <CoMISo/Utils/SmartPointer.hh>
#include <vector>

#include <adolc/adolc.h>
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <adolc/sparse/sparsedrivers.h>
#include <adolc/taping.h>

#include <Eigen/Eigen>
#if !(EIGEN_VERSION_AT_LEAST(3,1,0))
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#endif
#include <Eigen/Sparse>

#include <CoMISo/Config/CoMISoDefines.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

/** \class NProblemmInterfaceAD NProblemInterfaceAD.hpp <ACG/.../NProblemInterfaceAD.hh>

 Brief Description.

 Extend the problem interface with auto differentiation using ADOL-C
 */
class COMISODLLEXPORT NProblemInterfaceAD : public NProblemInterface
{
public:

    // Sparse Matrix Type
#if EIGEN_VERSION_AT_LEAST(3,1,0)
    typedef Eigen::SparseMatrix<double,Eigen::ColMajor> SMatrixNP;
#else
    typedef Eigen::DynamicSparseMatrix<double,Eigen::ColMajor> SMatrixNP;
#endif

    /// Default constructor
    NProblemInterfaceAD(int _n_unknowns) :
    n_unknowns_(_n_unknowns),
    x_d_(new adouble[_n_unknowns]),
    function_evaluated_(false),
    use_tape_(true),
    constant_hessian_evaluated_(false) {

        for(size_t i = 0; i < 11; ++i) tape_stats_[i] = -1;
    }

    /// Destructor
    virtual ~NProblemInterfaceAD() {
        delete[] x_d_;
    }

    // ================================================
    // Override these three methods
    // ================================================

    virtual void initial_x(double* _x) = 0;

    virtual adouble evaluate(const adouble* _x) = 0;

    virtual void store_result(const double* _x) = 0;

    // ================================================
    // Optionally override these methods, too
    // ================================================

    /**
     * \brief If the hessian is constant it only
     * has to be computed once if the function
     * is continuous over the entire domain and
     * taping is activated.
     */
    virtual bool constant_hessian() {
        return false;
    }

    /**
     * \brief Indicate whether the hessian is sparse.
     * If so, the computations (as well as the memory
     * consumption) can be performed more efficiently.
     */
    virtual bool sparse_hessian() {
        return false;
    }

    // ================================================

    virtual int n_unknowns() {

        return n_unknowns_;
    }

    virtual double eval_f(const double* _x) {

        adouble y_d = 0.0;
        double y = 0.0;

        trace_on(1); // Start taping

        // Fill data vector
        for(int i = 0; i < n_unknowns_; ++i) {
            x_d_[i] <<= _x[i];
        }

        // Call virtual function to compute
        // functional value
        y_d = evaluate(x_d_);

        y_d >>= y;

        trace_off();

#ifndef NDEBUG
        tapestats(1, tape_stats_);
        // Do some status output here...
#endif

        return y;
    }

    virtual void eval_gradient(const double* _x, double* _g) {

        if(!function_evaluated_ || !use_tape_) {
            // Evaluate original functional
            eval_f(_x);
            function_evaluated_ = true;
        }

        gradient(1, n_unknowns_, _x, _g);
    }

    virtual void eval_hessian(const double* _x, SMatrixNP& _H) {

        _H.resize(n_unknowns_, n_unknowns_);
        _H.setZero();

        if(constant_hessian() && constant_hessian_evaluated_) {
            _H = constant_hessian_;
            return;
        }

        if(!function_evaluated_ || !use_tape_) {
            // Evaluate original functional
            eval_f(_x);
            function_evaluated_ = true;
        }

        if(sparse_hessian()) {

            int nz = 0;
            int opt[2] = {0, 0};

            unsigned int* r_ind = NULL;
            unsigned int* c_ind = NULL;

            double* val = NULL;

            int ec = sparse_hess(1, n_unknowns_, 0, _x, &nz, &r_ind, &c_ind, &val, opt);

            assert(*nz >= 0);
            assert(r_ind != NULL);
            assert(c_ind != NULL);
            assert(val != NULL);

#ifndef NDEBUG
            std::cout << "Info: sparse_hessian() returned code " << ec << std::endl;
#endif

            for(int i = 0; i < nz; ++i) {

                _H.coeffRef(r_ind[i], c_ind[i]) = val[i];
            }

            if(constant_hessian()) {
                constant_hessian_ = _H;
                constant_hessian_evaluated_ = true;
            }

            delete[] val;
            delete[] r_ind;
            delete[] c_ind;

        } else {

            double** H;

            H = new double*[n_unknowns_];
            for(int i = 0; i < n_unknowns_; ++i) {
                H[i] = new double[i+1];
            }

            int ec = hessian(1, n_unknowns_, const_cast<double*>(_x), H);

#ifndef NDEBUG
            std::cout << "Info: hessian() returned code " << ec << std::endl;
#endif

            for(int i = 0; i < n_unknowns_; ++i) {
                for(int j = 0; j <= i; ++j) {

                    _H.coeffRef(i, j) = H[i][j];

                    if(i != j) {
                        _H.coeffRef(j, i) = H[i][j];
                    }
                }
            }

            if(constant_hessian()) {
                constant_hessian_ = _H;
                constant_hessian_evaluated_ = true;
            }

            for(int i = 0; i < n_unknowns_; ++i) {
                delete[] H[i];
            }
            delete[] H;
        }
    }

    virtual bool use_tape() const {
        return use_tape_;
    }

    /** \brief Use tape
     * Set this to false if the energy functional
     * is discontinuous (so that the operator tree
     * has to be re-established at each evaluation)
     */
    virtual void use_tape(bool _b) {
        use_tape_ = _b;
    }

private:

    int n_unknowns_;

    adouble* x_d_;

    int tape_stats_[11];

    bool function_evaluated_;
    bool use_tape_;

    SMatrixNP constant_hessian_;
    bool constant_hessian_evaluated_;
};

//=============================================================================
} // namespace COMISO
//=============================================================================
#endif // COMISO_ADOLC_AVAILABLE
//=============================================================================
#endif // COMISO_Eigen3_AVAILABLE
//=============================================================================
#endif /* NPROBLEMINTERFACEAD_HPP_ */
