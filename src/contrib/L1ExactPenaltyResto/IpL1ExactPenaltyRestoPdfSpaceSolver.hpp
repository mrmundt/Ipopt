//
// Created by David on 12/29/2020.
//

#ifndef SRC_IPL1EXACTPENALTYRESTOPDFSPACESOLVER_HPP
#define SRC_IPL1EXACTPENALTYRESTOPDFSPACESOLVER_HPP

#include "IpPDSystemSolver.hpp"
#include "IpAugSystemSolver.hpp"
#include "IpPDPerturbationHandler.hpp"

namespace Ipopt {
class L1ExactPenaltyRestoPDFSpaceSolver : public PDSystemSolver
{
public:
    L1ExactPenaltyRestoPDFSpaceSolver(
            AugSystemSolver& augSysSolver,
            PDPerturbationHandler& perturbHandler
    );
    ~L1ExactPenaltyRestoPDFSpaceSolver() override;

    bool InitializeImpl(
            const OptionsList& options,
            const std::string& prefix
            ) override;


    bool Solve(
            Number alpha,
            Number beta,
            const IteratesVector& rhs,
            IteratesVector& res,
            bool allow_inexact = false,
            bool improve_solution = false
            ) override;

    static void RegisterOptions(
            const SmartPtr<RegisteredOptions>& roptions
            );

private:
    L1ExactPenaltyRestoPDFSpaceSolver();

    L1ExactPenaltyRestoPDFSpaceSolver& operator=(
            const L1ExactPenaltyRestoPDFSpaceSolver&
            );

    SmartPtr<AugSystemSolver> augSysSolver_;
    SmartPtr<PDPerturbationHandler> perturbHandler_;
    CachedResults<void*> dummy_cache_;
    bool augsys_improved_{false};
    Index min_refinement_steps_{1};
    Index max_refinement_steps_{10}; // Defaults: from the option definitions
    Number residual_ratio_max_;
    Number residual_ratio_singular_;
    Number residual_improvement_factor_;
    Number neg_curv_test_tol_;
    bool neg_curv_test_reg_;


    bool SolveOnce(
            bool                  resolve_unmodified,
            bool                  pretend_singular,
            const SymMatrix&      W,
            const Matrix&         J_c,
            const Matrix&         J_d,
            const Matrix&         Px_L,
            const Matrix&         Px_U,
            const Matrix&         Pd_L,
            const Matrix&         Pd_U,
            const Vector&         z_L,
            const Vector&         z_U,
            const Vector&         v_L,
            const Vector&         v_U,
            const Vector&         slack_x_L,
            const Vector&         slack_x_U,
            const Vector&         slack_s_L,
            const Vector&         slack_s_U,
            const Vector&         sigma_x,
            const Vector&         sigma_s,
            Number                alpha,
            Number                beta,
            const IteratesVector& rhs,
            IteratesVector&       res
    );


    void ComputeResiduals(
            const SymMatrix&      W,
            const Matrix&         J_c,
            const Matrix&         J_d,
            const Matrix&         Px_L,
            const Matrix&         Px_U,
            const Matrix&         Pd_L,
            const Matrix&         Pd_U,
            const Vector&         z_L,
            const Vector&         z_U,
            const Vector&         v_L,
            const Vector&         v_U,
            const Vector&         slack_x_L,
            const Vector&         slack_x_U,
            const Vector&         slack_s_L,
            const Vector&         slack_s_U,
            const Vector&         sigma_x,
            const Vector&         sigma_s,
            Number                alpha,
            Number                beta,
            const IteratesVector& rhs,
            const IteratesVector& res,
            IteratesVector&       resid
    );


    Number ComputeResidualRatio(
            const IteratesVector& rhs,
            const IteratesVector& res,
            const IteratesVector& resid
    );

    void ScaleInvRhoRes(IteratesVector& res);
};
}

#endif //SRC_IPL1EXACTPENALTYRESTOPDFSPACESOLVER_HPP
