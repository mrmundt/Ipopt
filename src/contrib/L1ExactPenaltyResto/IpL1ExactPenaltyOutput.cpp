//
// Created by dav0 on 2/12/21.
//

#include "IpL1ExactPenaltyOutput.hpp"
#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"
#include "IpL1ExactPenaltyRestoData.hpp"

#include <cmath>


namespace Ipopt
{
L1ExactPenaltyOutput::L1ExactPenaltyOutput(
    const SmartPtr<OrigIterationOutput> &orig_iteration_output
)
        : orig_iteration_output_(orig_iteration_output)
        { }

L1ExactPenaltyOutput::~L1ExactPenaltyOutput() = default;

bool L1ExactPenaltyOutput::InitializeImpl(
        const OptionsList &options,
        const std::string &prefix)
{
    //options.GetBoolValue("print_info_string", print_info_string_, prefix);
    Index enum_int;
    options.GetEnumValue("inf_pr_output", enum_int, prefix);
    inf_pr_output_ = InfPrOutput(enum_int);
    options.GetIntegerValue("print_frequency_iter", print_frequency_iter_, prefix);
    options.GetNumericValue("print_frequency_time", print_frequency_time_, prefix);

    bool retval = true;
    if( IsValid(orig_iteration_output_))
    {
        retval = orig_iteration_output_->Initialize(Jnlst(),
                                                       IpNLP(),
                                                       IpData(),
                                                       IpCq(),
                                                       options,
                                                       prefix);
    }
    return retval;
}

void L1ExactPenaltyOutput::WriteOutput()
{
    auto l1epr_ipopt_nlp = dynamic_cast<const L1ExactPenaltyRestoIpoptNLP*>(&IpNLP());
    DBG_ASSERT(l1epr_ipopt_nlp);

    SmartPtr<IpoptData> orig_ip_data = &l1epr_ipopt_nlp->OrigIpData();
    SmartPtr<IpoptNLP> orig_ip_nlp = &l1epr_ipopt_nlp->OrigIpNLP();
    SmartPtr<IpoptCalculatedQuantities> orig_ip_cq = &l1epr_ipopt_nlp->OrigIpCq();

    Index iter = IpData().iter_count();
    orig_ip_data->Set_iter_count(iter);

    if( IsValid(orig_iteration_output_))
    {
        orig_iteration_output_->WriteOutput();
    }

    //std::string header = "iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls\n";
    std::string header = "iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_pr  ls rho\n";
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "*** Summary of Iteration %d for original NLP:", IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n\n");
    //if( IpData().info_iters_since_header() >= 10 && !IsValid(orig_iteration_output_))
    if( IpData().info_iters_since_header() >= 10)
    {
        Jnlst().Printf(J_ITERSUMMARY, J_MAIN,
                       header.c_str());
        IpData().Set_info_iters_since_header(0);
    }
    else
    {
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       header.c_str());
    }

    Number  inf_du = IpCq().curr_dual_infeasibility(NORM_MAX);

    Number mu = IpData().curr_mu();
    Number  dnrm = 0.;
    if(IsValid(IpData().delta()) && IsValid(IpData().delta()->x()) && IsValid(IpData().delta()->s()))
    {
        dnrm = Max(IpData().delta()->x()->Amax(), IpData().delta()->s()->Amax());
    }

    SmartPtr<const Vector> x = IpData().curr()->x();
    auto cx = dynamic_cast<const CompoundVector*>(GetRawPtr(x));
    DBG_ASSERT(cx);

    SmartPtr<const Vector> s = IpData().curr()->s();
    auto cs = dynamic_cast<const CompoundVector*>(GetRawPtr(s));
    DBG_ASSERT(cs);

    SmartPtr<IteratesVector> trial = orig_ip_data->trial()->MakeNewContainer();
    trial->Set_x(*cx->GetComp(0));
    trial->Set_s(*cs->GetComp(0));
    orig_ip_data->set_trial(trial);

    Number inf_pr = .0;
    switch (inf_pr_output_)
    {
        case IterationOutput::INTERNAL:
            inf_pr = orig_ip_cq->trial_primal_infeasibility(NORM_MAX);
            break;
        case IterationOutput::ORIGINAL:
            inf_pr = orig_ip_cq->unscaled_trial_nlp_constraint_violation(NORM_MAX);
            break;
    }
    Number f = orig_ip_cq->unscaled_trial_f();
    char info_iter = 'L';

    Number alpha_primal = IpData().info_alpha_primal();
    char alpha_primal_char = IpData().info_alpha_primal_char();
    Number alpha_dual = IpData().info_alpha_primal();
    Number regu_x = IpData().info_regu_x();
    char regu_x_buffer[8];
    char dashes[] = "   - ";
    char* regu_x_ptr;
    if( regu_x == .0)
    {
        regu_x_ptr = dashes;
    }
    else
    {
        Snprintf(regu_x_buffer, 7, "%5.1f", log10(regu_x));
        regu_x_ptr = regu_x_buffer;
    }
    Index ls_count = IpData().info_ls_count();
    const std::string info_string = IpData().info_string();

    Number current_time = 0.;
    Number last_output = IpData().info_last_output();
    // L1 Info

    auto l1_data = dynamic_cast<L1ExactPenaltyRestoData*>(&IpData().AdditionalData());
    DBG_ASSERT(l1_data);
    Number curr_rho = l1_data->GetCurrentRho();

    if((iter % print_frequency_iter_) == 0 && (print_frequency_time_ == 0. || last_output < (current_time = WallclockTime()) - print_frequency_time_ || last_output < .0))
    {
        //Jnlst().Printf(J_ITERSUMMARY, J_MAIN,
        //               "%4d%c%14.7e %7.2e %7.2e %5.1f %7.2e %5s %7.2e %7.2e%c%3d", iter, info_iter, f, inf_pr, inf_du, log10(mu), dnrm, regu_x_ptr, alpha_dual, alpha_primal, alpha_primal_char, ls_count);
        Jnlst().Printf(J_ITERSUMMARY, J_MAIN,
                       "%4d%c%14.7e %7.2e %7.2e %5.1f %7.2e %5s %7.2e%c%3d %7.2e", iter, info_iter, f, inf_pr, inf_du, log10(mu), dnrm, regu_x_ptr, alpha_primal, alpha_primal_char, ls_count, curr_rho);
        if(print_info_string_)
        {
            Jnlst().Printf(J_ITERSUMMARY, J_MAIN,
                           " %s", info_string.c_str());
        }
        else
        {
            Jnlst().Printf(J_DETAILED, J_MAIN,
                           " %s", info_string.c_str());
        }
        Jnlst().Printf(J_ITERSUMMARY, J_MAIN,
                       "\n");

        IpData().Set_info_last_output(current_time);
        IpData().Inc_info_iters_since_header();
    }

    if( Jnlst().ProduceOutput(J_DETAILED, J_MAIN))
    {
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "\n**************************************************\n");
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "*** Beginning Iteration %d from the following point:", IpData().iter_count());
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "\n**************************************************\n\n");

        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Primal infeasibility for restoration phase problem = %.16e\n", IpCq().curr_primal_infeasibility(NORM_MAX));
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Dual infeasibility for restoration phase problem   = %.16e\n", IpCq().curr_dual_infeasibility(NORM_MAX));

        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "||curr_x||_inf   = %.16e\n", IpData().curr()->x()->Amax());
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "||curr_s||_inf   = %.16e\n", IpData().curr()->s()->Amax());
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "||curr_y_c||_inf = %.16e\n", IpData().curr()->y_c()->Amax());
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "||curr_y_d||_inf = %.16e\n", IpData().curr()->y_d()->Amax());
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "||curr_z_L||_inf = %.16e\n", IpData().curr()->z_L()->Amax());
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "||curr_z_U||_inf = %.16e\n", IpData().curr()->z_U()->Amax());
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "||curr_v_L||_inf = %.16e\n", IpData().curr()->v_L()->Amax());
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "||curr_v_U||_inf = %.16e\n", IpData().curr()->v_U()->Amax());
    }
    if( Jnlst().ProduceOutput(J_MOREVECTOR, J_MAIN))
    {
        IpCq().curr_grad_lag_x()->Print(Jnlst(), J_MOREVECTOR, J_MAIN, "curr_grad_lag_x");
        IpCq().curr_grad_lag_s()->Print(Jnlst(), J_MOREVECTOR, J_MAIN, "curr_grad_lag_s");
        if(IsValid(IpData().delta()))
        {
            IpData().delta()->Print(Jnlst(), J_MOREVECTOR, J_MAIN, "delta");
        }
    }
    if(Jnlst().ProduceOutput(J_DETAILED, J_MAIN))
    {
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "\n\n***Current NLP Values for Iteration (L1 phase problem) %d:\n", IpData().iter_count());
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "\n                                   (scaled)                 (unscaled)\n");
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Objective...............: %24.16e  %24.16e\n", IpCq().curr_f(), IpCq().unscaled_curr_f());
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Dual infeasibility......: %24.16e  %24.16e\n", IpCq().curr_dual_infeasibility(NORM_MAX), IpCq().unscaled_curr_dual_infeasibility(NORM_MAX));
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Constraint violation....: %24.16e  %24.16e\n", IpCq().curr_nlp_constraint_violation(NORM_MAX), IpCq().unscaled_curr_nlp_constraint_violation(NORM_MAX));
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Complementarity.........: %24.16e  %24.16e\n", IpCq().curr_complementarity(0., NORM_MAX), IpCq().unscaled_curr_complementarity(0., NORM_MAX));
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Overall NLP error.......: %24.16e  %24.16e\n\n", IpCq().curr_nlp_error(), IpCq().unscaled_curr_nlp_error());
    }
    if( Jnlst().ProduceOutput(J_VECTOR, J_MAIN) )
    {
        IpCq().curr_grad_f()->Print(Jnlst(), J_VECTOR, J_MAIN, "grad_f");
        IpCq().curr_c()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_c");
        IpCq().curr_d()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_d");
        IpCq().curr_d_minus_s()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_d - curr_s");
    }

    if( Jnlst().ProduceOutput(J_MATRIX, J_MAIN) )
    {
        IpCq().curr_jac_c()->Print(Jnlst(), J_MATRIX, J_MAIN, "jac_c");
        IpCq().curr_jac_d()->Print(Jnlst(), J_MATRIX, J_MAIN, "jac_d");
        IpData().W()->Print(Jnlst(), J_MATRIX, J_MAIN, "W");
    }

    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n\n");
}


}