//=============================================================================
//
//  CLASS NConstraintInterfaceAD
//
//=============================================================================

#ifndef COMISO_NCONSTRAINTINTERFACEAD_HH
#define COMISO_NCONSTRAINTINTERFACEAD_HH

//== COMPILE-TIME PACKAGE REQUIREMENTS ========================================
#include <CoMISo/Config/config.hh>
#if COMISO_ADOLC_AVAILABLE
#if COMISO_Eigen3_AVAILABLE

//== INCLUDES =================================================================

#include <CoMISo/Config/CoMISoDefines.hh>
#include "SuperSparseMatrixT.hh"

#include <adolc/adolc.h>
#include <adolc/adouble.h>
#include <adolc/drivers/drivers.h>
#include <adolc/sparse/sparsedrivers.h>
#include <adolc/taping.h>

#include "NConstraintInterface.hh"

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace COMISO {

//== CLASS DEFINITION =========================================================

/** \class NProblemInterfaceAD NProblemInterfaceAD.hpp

    The problem interface using automatic differentiation.
 */
class COMISODLLEXPORT NConstraintInterfaceAD : public NConstraintInterface {
public:

    // Define Sparse Datatypes
    typedef NConstraintInterface::SVectorNC SVectorNC;
    typedef NConstraintInterface::SMatrixNC SMatrixNC;

    typedef NConstraintInterface::ConstraintType ConstraintType;

    /// Default constructor
    NConstraintInterfaceAD(int _n_unknowns, const ConstraintType _type = NC_EQUAL) :
        NConstraintInterface(_type),
        n_unknowns_(_n_unknowns),
        type_(_type),
        x_d_(new adouble[n_unknowns_]),
        r_ind_(NULL),
        c_ind_(NULL),
        val_(NULL),
        H_(NULL),
        grad_(new double[n_unknowns_]),
        function_evaluated_(false),
        use_tape_(true),
        constant_hessian_evaluated_(false) {

        for(size_t i = 0; i < 11; ++i) tape_stats_[i] = -1;

        if(sparse_hessian()) {
            r_ind_ = new unsigned int[n_unknowns_];
            c_ind_ = new unsigned int[n_unknowns_];
            val_ = new double[n_unknowns_];
        } else {
            H_ = new double*[n_unknowns_];
            for(int i = 0; i < n_unknowns_; ++i) {
                H_[i] = new double[i+1];
            }
        }
    }

    /// Destructor
    virtual ~NConstraintInterfaceAD() {
        delete[] grad_;
        delete[] x_d_;

        if(sparse_hessian()) {
            delete[] r_ind_;
            delete[] c_ind_;
            delete[] val_;
        } else {
            for(int i = 0; i < n_unknowns_; ++i) {
                delete[] H_[i];
            }
            delete[] H_;
        }
    }

    /**
     * \brief Only override this function
     */
    virtual adouble evaluate(const adouble* _x) = 0;

public:

    virtual int n_unknowns() {
        return n_unknowns_;
    }

    virtual double eval_constraint(const double* _x) {

        double y = 0.0;

        if(!function_evaluated_ || !use_tape_) {

            adouble y_d = 0.0;

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

            function_evaluated_ = true;

        } else {

            double ay[1] = {0.0};

            int ec = function(1, 1, n_unknowns_, const_cast<double*>(_x), ay);

#ifndef NDEBUG
            std::cout << "Info: function() returned code " << ec << std::endl;
#endif

            y = ay[0];
        }

        return y;
    }

    virtual void eval_gradient(const double* _x, SVectorNC& _g) {

        if(!function_evaluated_ || !use_tape_) {
            // Evaluate original functional
            eval_constraint(_x);
        }

        _g.resize(n_unknowns_);
        _g.setZero();

        int ec = gradient(1, n_unknowns_, _x, grad_);

#ifndef NDEBUG
        std::cout << "Info: gradient() returned code " << ec << std::endl;
#endif

        for(int i = 0; i < n_unknowns_; ++i) {
            _g.coeffRef(i) = grad_[i];
        }
    }

    virtual void eval_hessian(const double* _x, SMatrixNC& _H) {

        _H.resize(n_unknowns_, n_unknowns_);

        if(constant_hessian() && constant_hessian_evaluated_) {
            _H = constant_hessian_;
            return;
        }

        if(!function_evaluated_ || !use_tape_) {
            // Evaluate original functional
            eval_constraint(_x);
        }

        if(sparse_hessian()) {

            int nz = 0;
            int opt[2] = {0, 0};

            int ec = sparse_hess(1, n_unknowns_, 0, _x, &nz, &r_ind_, &c_ind_, &val_, opt);

            assert(*nz >= 0);
            assert(r_ind != NULL);
            assert(c_ind != NULL);
            assert(val != NULL);

#ifndef NDEBUG
            std::cout << "Info: sparse_hessian() returned code " << ec << std::endl;
#endif

            for(int i = 0; i < nz; ++i) {

                _H(r_ind_[i], c_ind_[i]) = val_[i];
            }

            if(constant_hessian()) {
                constant_hessian_ = _H;
                constant_hessian_evaluated_ = true;
            }

        } else {

            int ec = hessian(1, n_unknowns_, const_cast<double*>(_x), H_);

#ifndef NDEBUG
            std::cout << "Info: hessian() returned code " << ec << std::endl;
#endif

            for(int i = 0; i < n_unknowns_; ++i) {
                for(int j = 0; j <= i; ++j) {

                    _H(i, j) = H_[i][j];

                    if(i != j) {
                        _H(j, i) = H_[i][j];
                    }
                }
            }

            if(constant_hessian()) {
                constant_hessian_ = _H;
                constant_hessian_evaluated_ = true;
            }
        }
    }

    /** \brief Use tape
     * Set this to false if the energy functional
     * is discontinuous (so that the operator tree
     * has to be re-established at each evaluation)
     */
    virtual bool use_tape() const {
        return use_tape_;
    }

    virtual void use_tape(bool _b) {
        use_tape_ = _b;
    }

    /**
     * \brief Override this function and let it return
     * true if you want to speed up the computation
     * by exploiting the fact that the hessian
     * matrix is constant and thus has to be evaluated
     * only once.
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

private:

    // Number of unknowns
    int n_unknowns_;

    // Constraint type
    ConstraintType type_;

    adouble* x_d_;

    // Sparse hessian data
    unsigned int* r_ind_;
    unsigned int* c_ind_;
    double* val_;

    // Non-sparse hessian data
    double** H_;

    double* grad_;

    int tape_stats_[11];

    bool function_evaluated_;
    bool use_tape_;

    SMatrixNC constant_hessian_;
    bool constant_hessian_evaluated_;
};

//=============================================================================
}// namespace COMISO
//=============================================================================
#endif // COMISO_ADOLC_AVAILABLE
//=============================================================================
#endif // COMISO_Eigen3_AVAILABLE
//=============================================================================
#endif // ACG_NCONSTRAINTINTERFACEAD_HH defined
//=============================================================================

