//
// Created by David on 1/6/2021.
//

#ifndef SRC_IPL1IPOPTALG_HPP
#define SRC_IPL1IPOPTALG_HPP

#include "IpIpoptNLP.hpp"
#include "IpAlgStrategy.hpp"
#include "IpSearchDirCalculator.hpp"
#include "IpLineSearch.hpp"
#include "IpMuUpdate.hpp"
#include "IpConvCheck.hpp"
#include "IpOptionsList.hpp"
#include "IpIterateInitializer.hpp"
#include "IpIterationOutput.hpp"
#include "IpAlgTypes.hpp"
#include "IpHessianUpdater.hpp"
#include "IpEqMultCalculator.hpp"
#include "IpL1ExactPenaltyRhoUpdater.hpp"

namespace Ipopt
{
DECLARE_STD_EXCEPTION(L1_STEP_COMPUTATION_FAILED);

class L1IpoptAlg : public AlgorithmStrategyObject
{
public:
    L1IpoptAlg(
            const SmartPtr<SearchDirectionCalculator>& search_dir_calculator,
            const SmartPtr<LineSearch>&                line_search,
            const SmartPtr<MuUpdate>&                  mu_update,
            const SmartPtr<ConvergenceCheck>&          conv_check,
            const SmartPtr<IterateInitializer>&        iterate_initializer,
            const SmartPtr<IterationOutput>&           iter_output,
            const SmartPtr<HessianUpdater>&            hessian_updater,
            const SmartPtr<L1ExactPenaltyRhoUpdater>&  l1exactpenalty_rho_updater,
            const SmartPtr<EqMultiplierCalculator>&    eq_multiplier_calculator = NULL);

    ~L1IpoptAlg() override;

    bool InitializeImpl(
            const OptionsList& options,
            const std::string& prefix
            ) override;

    SolverReturn Optimize(
            bool isResto = false
            );

    static void RegisterOptions(
            SmartPtr<RegisteredOptions> roptions
            );

    static void print_copyright_message(
            const Journalist& jnlst
            );

private:
    L1IpoptAlg();

    L1IpoptAlg(
            const L1IpoptAlg&
            );

    SmartPtr<SearchDirectionCalculator> search_dir_calculator_;
    SmartPtr<LineSearch> line_search_;
    SmartPtr<MuUpdate> mu_update_;
    SmartPtr<ConvergenceCheck> conv_check_;
    SmartPtr<IterateInitializer> iterate_initializer_;
    SmartPtr<IterationOutput> iter_ouput_;
    SmartPtr<HessianUpdater> hessian_updater_;

    SmartPtr<EqMultiplierCalculator> eq_multiplier_calculator_;
    SmartPtr<L1ExactPenaltyRhoUpdater> l1exactpenalty_rho_updater_;

    void UpdateHessian();
    bool UpdateBarrierParameter();
    bool ComputeSearchDirection();
    void ComputeAcceptableTrialPoint();
    void AcceptTrialPoint();
    void OutputIteration();
    void InitializeIterates();
    void PrintProblemStatistics();
    void ComputeFeasibilityMultipliers();
    bool skip_print_problem_stats_;
    Number kappa_sigma_;
    bool recalc_y_;
    Number recalc_y_feas_tol_;
    bool mehrotra_algorithm_;
    std::string linear_solver_;

    void calc_number_of_bounds(
            const Vector& x,
            const Vector& x_L,
            const Vector& x_U,
            const Matrix& Px_L,
            const Matrix& Px_U,
            Index&        n_tot,
            Index&        n_only_lower,
            Index&        n_both,
            Index&        n_only_upper
    );

    Number correct_bound_multiplier(
            const Vector&           trial_z,
            const Vector&           trial_slack,
            const Vector&           trial_compl,
            SmartPtr<const Vector>& new_trial_z
    );

    void ComputeRhoTrial();
    void UpdateRhoAction();
};
}

#endif //SRC_IPL1IPOPTALG_HPP
