//
// Created by David on 1/3/2021.
//

#include "IpL1ExactPenaltyRestoFilterConvCheck.hpp"
#include "IpCompoundVector.hpp"
//#include "IpRestoIpoptNLP.hpp"
#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"
#include "IpRestoPhase.hpp"

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
    static const Index dbg_verbosity = 0;
#endif

L1ExactPenaltyRestoFilterConvCheck::L1ExactPenaltyRestoFilterConvCheck()
{
    DBG_START_FUN("L1ExactPenaltyRestoFilterConvCheck::L1ExactPenaltyRestoFilterConvCheck()",
                  dbg_verbosity);
}

L1ExactPenaltyRestoFilterConvCheck::~L1ExactPenaltyRestoFilterConvCheck() noexcept
{
    DBG_START_FUN("L1ExactPenaltyRestoFilterConvCheck::~L1ExactPenaltyRestoFilterConvCheck()",
                  dbg_verbosity);
}

void L1ExactPenaltyRestoFilterConvCheck::RegisterOptions(
        SmartPtr<RegisteredOptions> roptions)
{}

bool L1ExactPenaltyRestoFilterConvCheck::InitializeImpl(
        const OptionsList &options, const std::string &prefix)
{
    return OptimalityErrorConvergenceCheck::InitializeImpl(options, prefix);
}

ConvergenceCheck::ConvergenceStatus L1ExactPenaltyRestoFilterConvCheck::CheckConvergence(
        bool call_intermediate_callback)
{
    const L1ExactPenaltyRestoIpoptNLP* l1epr_ipopt_nlp = static_cast<const L1ExactPenaltyRestoIpoptNLP*>(&IpNLP());
    DBG_ASSERT(dynamic_cast<const L1ExactPenaltyRestoIpoptNLP*>(&IpNLP()))

    SmartPtr<IpoptData> orig_ip_data = &l1epr_ipopt_nlp->OrigIpData();
    SmartPtr<IpoptCalculatedQuantities> orig_ip_cq = &l1epr_ipopt_nlp->OrigIpCq();

    SmartPtr<const Vector> x = IpData().curr()->x();
    const CompoundVector* cx = static_cast<const CompoundVector*>(GetRawPtr(x));
    DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(x)));
    SmartPtr<IteratesVector> trial = orig_ip_data->curr()->MakeNewContainer();
    trial->Set_x(*cx->GetComp(0));
    trial->Set_s(IpData().curr()->s());
    orig_ip_data->set_trial(trial);

    if(call_intermediate_callback) {
        AlgorithmMode mode = RestorationPhaseMode;

        Index iter = IpData().iter_count();
        Number inf_pr = orig_ip_cq->trial_primal_infeasibility(NORM_MAX);
        Number inf_du = IpCq().curr_dual_infeasibility(NORM_MAX);
        Number mu = IpData().curr_mu();
        Number dnrm;
        if (IsValid(IpData().delta()) && IsValid(IpData().delta()->x()) &&
            IsValid(IpData().delta()->s()))
        {
            dnrm = Max(IpData().delta()->x()->Amax(), IpData().delta()->s()->Amax());
        }
        else
        {
            dnrm = 0.0;
        }
        Number alpha_primal = IpData().info_alpha_primal();
        Number alpha_dual = IpData().info_alpha_dual();
        Number regu_x = IpData().info_regu_x();
        Number unscaled_f = orig_ip_cq->unscaled_trial_f();
        Index ls_count = IpData().info_ls_count();
        bool request_stop = !IpNLP().IntermediateCallBack(mode, iter, unscaled_f, inf_pr, inf_du, du, dnrm, regu_x,
                                                          alpha_dual, alpha_primal, ls_count, &IpData(), &IpCq());

        if( request_stop )
        {
            return ConvergenceCheck::USER_STOP;
        }
    }

    if (IpData().iter_count() >= maximum_iters_)
    {
        return ConvergenceCheck::MAXITER_EXCEEDED;
    }

    if( succesive_resto_iter_ > maximum_resto_iters_)
    {
        Jnlst().Printf(J_WARNING, J_MAIN,
                       "Max resto iterations is skipped.")
    }
    succesive_resto_iter_++;

    // We do not check the original filter.

    ConvergenceStatus status = CONTINUE;

    // Calculate the f and theta for the original problem
    Number orig_trial_theta = orig_ip_cq->trial_constraint_violation();
    Number orig_curr_theta = orig_ip_cq->curr_constraint_violation();

    // check acceptability to the filter
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "orig_curr_theta = %8.2e, orig_trial_theta = %8.2e\n", orig_curr_theta, orig_trial_theta);


    Number orig_curr_inf_pr = orig_ip_cq->curr_primal_infeasibility(NORM_MAX);
    Number orig_trial_inf_pr = orig_ip_cq->trial_primal_infeasibility(NORM_MAX);
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "orig_curr_inf_pr = %8.2e, orig_trial_inf_pr = %8.2e\n", orig_curr_inf_pr, orig_trial_inf_pr);

    Number orig_inf_pr_max = Max(kappa_resto_ * orig_curr_inf_pr, Min(orig_ip_data->tol(), orig_constr_viol_tol_));
    if( kappa_resto_ == 0. )
    {
        orig_inf_pr_max = 0.;
    }

    if( first_resto_iter_ )
    {
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "This is the first iteration - continue to take at least one step.\n");
        status = CONTINUE;
    }
    else if( orig_ip_cq->IsSquareProblem() && orig_trial_inf_pr <= Min(orig_ip_data->tol(), orig_constr_viol_tol_) )
    {
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Restoration phase found points satisfying feasibility tolerance in square problem.\n");
        status = CONVERGED;
    }
    //else if( orig_trial_inf_pr > orig_inf_pr_max )
    //{
        // One could have a minimum infeasibility solution otherwise.
        // I think one wishes to only check the original filter if we have
        // attained certain progress toward feasibility.
        //Jnlst().Printf(J_DETAILED, J_MAIN,
        //               "Point does not provide sufficient reduction w.r.t the original constraint violation (orig_inf_pr_max=%e).\n",
        //               orig_inf_pr_max);
        //status = CONTINUE;
        // I want it to continue anyways.
    //}

    if( status == CONTINUE )
    {
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Checking convergence for restoration phase problem...\n");
        status = OptimalityErrorConvergenceCheck::CheckConvergence(false);
        if( status == CONVERGED || status == CONVERGED_TO_ACCEPTABLE_POINT )
        {
            Number orig_trial_primal_inf = orig_ip_cq->trial_primal_infeasibility(NORM_MAX);
            // ToDo make the factor in following line an option
            if( orig_trial_primal_inf <= 1e2 * IpData().tol() )
            {
                //        if (orig_trial_primal_inf <= 1e2*orig_ip_data->tol()) {
                if( IpData().tol() > 1e-1 * orig_ip_data->tol() )
                {
                    // For once, we tighten the convergence tolerance for the
                    // restoration phase problem in case the problem is only
                    // very slightly infeasible.
                    IpData().Set_tol(1e-2 * IpData().tol());
                    status = CONTINUE;
                    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                                   "Tightening restoration phase tolerance to %e.\n", IpData().tol());
                    IpData().Append_info_string("!");
                }
                else
                {
                    Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                                   "Restoration phase converged to a feasible point that is\n"
                                   "unacceptable to the filter for the original problem.\n");
                    // I just want it to return "CONVERGED"
                    //THROW_EXCEPTION(RESTORATION_CONVERGED_TO_FEASIBLE_POINT,
                    //                "Restoration phase converged to a feasible point that is "
                    //                "unacceptable to the filter for the original problem.");
                }
            }
            else
            {
                THROW_EXCEPTION(LOCALLY_INFEASIBLE, "Restoration phase converged to a point of local infeasibility");
            }
        }
    }

}


}