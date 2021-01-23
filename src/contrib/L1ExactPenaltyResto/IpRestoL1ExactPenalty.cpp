//
// Created by David on 1/14/2021.
//

#include "IpRestoL1ExactPenalty.hpp"
#include "IpL1ExactPenaltyRestoCalculatedQuantities.hpp"
#include "IpDefaultIterateInitializer.hpp"

namespace Ipopt
{
L1ExactPenaltyRestorationPhase::L1ExactPenaltyRestorationPhase(L1IpoptAlg &resto_alg,
                                                               const SmartPtr<EqMultiplierCalculator> &eq_mult_calculator) :
                                         resto_alg_(&resto_alg),
                                         eq_mult_calculator_(eq_mult_calculator)
{}

bool L1ExactPenaltyRestorationPhase::PerformRestoration()
{
    DBG_START_METH("L1ExactPenaltyRestorationPhase::PerformRestoration",
                   dbg_verbosity);
    DBG_ASSERT(IpCq().curr_constraint_violation() > 0.);
    SmartPtr<IpoptAdditionalData> l1data = new L1ExactPenaltyRestoData();
    SmartPtr<IpoptData> l1_ip_data = new IpoptData(l1data, IpData().cpu_time_start());
    SmartPtr<IpoptNLP> l1_ip_nlp = new L1ExactPenaltyRestoIpoptNLP(
            IpNLP(),
            IpData(),
            IpCq(),
            l1_ip_data);
    // The actual calculated quantities for the l1 mode.
    SmartPtr<IpoptCalculatedQuantities> l1_ip_cq = new L1ExactPenaltyRestoCQ(l1_ip_nlp, l1_ip_data);

    resto_alg_->Initialize(Jnlst(),
                           *l1_ip_nlp,
                           *l1_ip_data,
                           *l1_ip_cq,
                           *resto_options_,
                           "resto.");

    l1_ip_data->Set_iter_count(IpData().iter_count() + 1);
    l1_ip_data->Set_info_regu_x(IpData().info_regu_x());
    l1_ip_data->Set_info_alpha_primal(IpData().info_alpha_primal());
    l1_ip_data->Set_info_alpha_primal_char(IpData().info_alpha_primal_char());
    l1_ip_data->Set_info_alpha_dual(IpData().info_alpha_dual());
    l1_ip_data->Set_info_ls_count(IpData().info_ls_count());
    l1_ip_data->Set_info_iters_since_header(IpData().info_iters_since_header());
    l1_ip_data->Set_info_last_output(IpData().info_last_output());

    SolverReturn l1_status = resto_alg_->Optimize(true);

    int retval = -1;

    if( l1_status != SUCCESS)
    {
        SmartPtr<const IteratesVector> l1_curr = l1_ip_data->curr();
        if(IsValid(l1_curr))
        {
            SmartPtr<IteratesVector> trial = IpData().trial()->MakeNewContainer();
            // I use dynamic casting, hopefully this won't be an issue.
            SmartPtr<const Vector> l1_curr_x = l1_curr->x();
            auto cx = dynamic_cast<const CompoundVector*>(GetRawPtr(l1_curr_x));
            DBG_ASSERT(cx);

            SmartPtr<const Vector> l1_curr_s = l1_curr->s();
            auto cs = dynamic_cast<const CompoundVector*>(GetRawPtr(l1_curr_s));

            SmartPtr<const Vector> l1_curr_y_c = l1_curr->y_c();
            auto cy_c = dynamic_cast<const CompoundVector*>(GetRawPtr(l1_curr_y_c));

            SmartPtr<const Vector> l1_curr_y_d = l1_curr->y_d();
            auto cy_d = dynamic_cast<const CompoundVector*>(GetRawPtr(l1_curr_y_d));

            SmartPtr<const Vector> l1_curr_z_L = l1_curr->z_L();
            auto c_zL = dynamic_cast<const CompoundVector*>(GetRawPtr(l1_curr_z_L));

            SmartPtr<const Vector> l1_curr_z_U = l1_curr->z_U();
            auto c_zU = dynamic_cast<const CompoundVector*>(GetRawPtr(l1_curr_z_U));

            SmartPtr<const Vector> l1_curr_v_L = l1_curr->v_L();
            auto c_vL = dynamic_cast<const CompoundVector*>(GetRawPtr(l1_curr_v_L));

            SmartPtr<const Vector> l1_curr_v_U = l1_curr->v_U();
            auto c_vU = dynamic_cast<const CompoundVector*>(GetRawPtr(l1_curr_v_U));

            trial->Set_primal(*cx->GetComp(0), *cs->GetComp(0));
            trial->Set_eq_mult(*cy_c->GetComp(0), *cy_d->GetComp(0));
            trial->Set_bound_mult(*c_zL->GetComp(0), *c_zU->GetComp(0), *c_vL->GetComp(0), *c_vU->GetComp(0));
            IpData().set_trial(trial);
            IpData().AcceptTrialPoint();
        }
    }

    if(l1_status == SUCCESS)
    {
        if(Jnlst().ProduceOutput(J_DETAILED, J_LINE_SEARCH))
        {
            Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                           "\nRESTORATION PHASE RESULTS\n");
            Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                           "\n\nOptimal solution found! \n");
            Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                           "Optimal Objective Value = %.16E\n", l1_ip_cq->curr_f());
            Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                           "Number of Iterations = %d\n", l1_ip_data->iter_count());
        }
        if(Jnlst().ProduceOutput(J_VECTOR, J_LINE_SEARCH))
        {
            l1_ip_data->curr()->Print(Jnlst(), J_VECTOR, J_LINE_SEARCH, "curr");
        }
        retval = 0;
    }
    else if( l1_status == STOP_AT_TINY_STEP || l1_status == STOP_AT_ACCEPTABLE_POINT)
    {
        Number orig_primal_inf = IpCq().curr_primal_infeasibility(NORM_MAX);

        if( orig_primal_inf <= IpData().tol() )
        {
            Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                           "Restoration phase converged to a point with small primal infeasibility.\n"
                           "Original primal inf is less than the tolerance.");
            //THROW_EXCEPTION(RESTORATION_CONVERGED_TO_FEASIBLE_POINT,
            //                "Restoration phase converged to a point with small primal infeasibility");
            retval = 0;
        }
        else
        {
            THROW_EXCEPTION(LOCALLY_INFEASIBLE, "Restoration phase converged to a point of local infeasibility");
        }
    }
    else if ( l1_status == MAXITER_EXCEEDED)
    {
        THROW_EXCEPTION(RESTORATION_MAXITER_EXCEEDED,
                        "Maximal number of iteration exceeded in restoration phase.");
    }
    else if (l1_status == CPUTIME_EXCEEDED)
    {
        THROW_EXCEPTION(RESTORATION_CPUTIME_EXCEEDED,
                        "Maximal CPU time exceeded in restoration phase.");
    }
    else if (l1_status == LOCAL_INFEASIBILITY)
    {
        THROW_EXCEPTION(LOCALLY_INFEASIBLE,
                        "Restoration phase converged to a point of local infeasibility");
    }
    else if (l1_status == RESTORATION_FAILURE)
    {
        Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                       "Restoration phase in the restoraion phase failed.");
        THROW_EXCEPTION(RESTORATION_FAILED, "Restoration phase in the restoration phase failed.")
    }
    else if (l1_status == USER_REQUESTED_STOP)
    {
        THROW_EXCEPTION(RESTORATION_USER_STOP, "User requested stop during restoration phase.");
    }
    else
    {
        Jnlst().Printf(J_ERROR, J_MAIN,
                       "Unknown return status.\n I.e. not SUCCESS, MAXITER_EXCEED, RESTO_CPU_EXCEED, LOCAL_INFES, RESTO_FAIL so far.");
        retval = 1;
    }
    if (retval == 0)
    {
        SmartPtr<const Vector> l1_curr_x = l1_ip_data->curr()->x();
        SmartPtr<const CompoundVector> cx = dynamic_cast<const CompoundVector*>(GetRawPtr(l1_curr_x));
        DBG_ASSERT(cx);

        SmartPtr<const Vector> l1_curr_s = l1_ip_data->curr()->s();
        SmartPtr<const CompoundVector> cs = dynamic_cast<const CompoundVector*>(GetRawPtr(l1_curr_s));

        SmartPtr<IteratesVector> trial = IpData().trial()->MakeNewContainer();
        trial->Set_primal(*cx->GetComp(0), *cs->GetComp(0));
        IpData().set_trial(trial);

        if (IpCq().IsSquareProblem())
        {
            Number const_viol = IpCq().unscaled_curr_nlp_constraint_violation(NORM_MAX);

            if(const_viol <= constr_viol_tol_)
            {
                Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                               "Recursive restoration phase algorithm terminated successfully for square problem.\n");
                IpData().AcceptTrialPoint();
                THROW_EXCEPTION(FEASIBILITY_PROBLEM_SOLVED,
                                "Restoration phase converged to sufficiently feasible point of original square problem\n")
            }
        }

        SmartPtr<IteratesVector> delta = IpData().curr()->MakeNewIteratesVector(true);
        delta->Set(0.);
        ComputeBoundMultiplierStep(*delta->z_L_NonConst(),
                                   *IpData().curr()->z_L(),
                                   *IpCq().curr_slack_x_L(),
                                   *IpCq().trial_slack_x_L());
        ComputeBoundMultiplierStep(*delta->z_U_NonConst(),
                                   *IpData().curr()->z_U(),
                                   *IpCq().curr_slack_x_U(),
                                   *IpCq().trial_slack_s_L());
        ComputeBoundMultiplierStep(*delta->v_L_NonConst(),
                                   *IpData().curr()->v_L(),
                                   *IpCq().curr_slack_s_L(),
                                   *IpCq().trial_slack_s_L());
        ComputeBoundMultiplierStep(*delta->v_U_NonConst(),
                                   *IpData().curr()->v_U(),
                                   *IpCq().curr_slack_s_U(),
                                   *IpCq().trial_slack_s_U());

        DBG_PRINT_VECTOR(1, "delta_z_L", *delta->z_L());
        DBG_PRINT_VECTOR(1, "delta_z_U", *delta->z_U());
        DBG_PRINT_VECTOR(1, "delta_v_L", *delta->v_L());
        DBG_PRINT_VECTOR(1, "delta_v_U", *delta->v_U());

        Number alpha_dual = IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
                                                          *delta->z_L_NonConst(),
                                                          *delta->z_U_NonConst(),
                                                          *delta->v_L_NonConst(),
                                                          *delta->v_U_NonConst());

        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "Step size for bound multipliers: %8.2\n", alpha_dual);

        IpData().SetTrialBoundMultipliersFromStep(alpha_dual,
                                                  *delta->z_L(),
                                                  *delta->z_U(),
                                                  *delta->v_L(),
                                                  *delta->v_U());

        Number bound_mult_max = Max(IpData().trial()->z_L()->Amax(),
                                    IpData().trial()->z_U()->Amax(),
                                    IpData().trial()->v_L()->Amax(),
                                    IpData().trial()->v_U()->Amax());
        if (bound_mult_max > bound_mult_reset_threshold_)
        {
            trial = IpData().trial()->MakeNewContainer();
            Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                           "Bound multiplier after restoration phase too large (max=%8.2e). Set all to 1.\n", bound_mult_max);
            trial->create_new_z_L();
            trial->create_new_z_U();
            trial->create_new_v_L();
            trial->create_new_v_U();
            trial->z_L_NonConst()->Set(1.);
            trial->z_U_NonConst()->Set(1.);
            trial->v_L_NonConst()->Set(1.);
            trial->v_U_NonConst()->Set(1.);
            IpData().set_trial(trial);
        }

        DefaultIterateInitializer::least_square_mults(Jnlst(),
                                                      IpNLP(),
                                                      IpData(),
                                                      IpCq(),
                                                      eq_mult_calculator_,
                                                      constr_mult_reset_threshold_);

        DBG_PRINT_VECTOR(2, "y_c", *IpData().curr()->y_c());
        DBG_PRINT_VECTOR(2, "y_d", *IpData().curr()->y_d());

        IpData().Set_iter_count(l1_ip_data->iter_count() - 1 );

        IpData().Set_info_skip_output(true);
        IpData().Set_info_iters_since_header(l1_ip_data->info_iters_since_header());
        IpData().Set_info_last_output(l1_ip_data->info_last_output());

    }
    return (retval == 0);
}

void L1ExactPenaltyRestorationPhase::ComputeBoundMultiplierStep(Vector &delta_z,
                                                                const Vector &curr_z,
                                                                const Vector &curr_slack,
                                                                const Vector &trial_slack)
{
    Number mu = IpData().curr_mu();

    delta_z.Copy(curr_slack);
    delta_z.Axpy(-1., trial_slack);
    delta_z.ElementWiseMultiply(curr_z);
    delta_z.AddScalar(mu);
    delta_z.ElementWiseDivide(curr_slack);
    delta_z.Axpy(-1., curr_z);
}


}
