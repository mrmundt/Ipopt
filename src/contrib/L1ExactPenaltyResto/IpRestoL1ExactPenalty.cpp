//
// Created by David on 1/14/2021.
//

#include "IpRestoL1ExactPenalty.hpp"
#include "IpL1ExactPenaltyRestoData.hpp"
#include "IpL1ExactPenaltyRestoCalculatedQuantities.hpp"

namespace Ipopt
{
RestoL1ExactPenalty::RestoL1ExactPenalty(IpL1IpoptAlg &resto_alg,
                                         const SmartPtr<EqMultiplierCalculator> &eq_mult_calculator) :
                                         resto_alg_(&resto_alg),
                                         eq_mult_calculator_(eq_mult_calculator)
{}

bool RestoL1ExactPenalty::PerformRestoration()
{
    DBG_START_METH("RestoL1ExactPenalty::PerformRestoration",
                   dbg_verbosity);
    DBG_ASSERT(IpCq().curr_constraint_violation() > 0.);
    SmartPtr<IpoptAdditionalData> l1data = new L1ExactPenaltyRestoData();
    SmartPtr<IpoptData> l1_ip_data = new IpoptData(l1data, IpData().cpu_time_start());
    SmartPtr<IpoptNLP> l1_ip_nlp = new L1ExactPenaltyRestoIpoptNLP(IpNLP(), IpData(), IpCq());
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
}
}
