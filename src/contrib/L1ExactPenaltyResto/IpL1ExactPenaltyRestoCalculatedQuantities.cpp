//
// Created by David on 10/12/2020.
//

#include "IpL1ExactPenaltyRestoCalculatedQuantities.hpp"

#include <cmath>
#include <limits>


namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
    static const Index dbg_verbosity = 0;
#endif

L1ExactPenaltyRestoCQ::L1ExactPenaltyRestoCQ(
        const SmartPtr<IpoptNLP>& ip_nlp,
        const SmartPtr<IpoptData>& ip_data
        )
        : IpoptCalculatedQuantities(ip_nlp, ip_data),
        ip_nlp_l1_(ip_nlp),
        ip_data_l1_(ip_data),
        curr_f_cache_l1_(1),
        dampind_x_L_cache_l1_(2),
        dampind_x_U_cache_l1_(2),
        dampind_s_L_cache_l1_(2),
        dampind_s_U_cache_l1_(2),
        trial_f_cache_l1_(5),
        curr_grad_f_cache_l1_(2),
        trial_grad_f_cache_l1_(1),
        curr_barrier_obj_cache_l1_(2),
        trial_barrier_obj_cache_l1_(5),
        curr_grad_barrier_obj_x_cache_l1_(1),
        curr_grad_barrier_obj_s_cache_l1_(1),
        grad_kappa_times_damping_x_cache_l1_(1),
        grad_kappa_times_damping_s_cache_l1_(1),
        curr_compl_x_L_cache_l1_(1),
        curr_compl_x_U_cache_l1_(1),
        curr_compl_s_L_cache_l1_(1),
        curr_compl_s_U_cache_l1_(1),
        trial_compl_x_L_cache_l1_(1),
        trial_compl_x_U_cache_l1_(1),
        trial_compl_s_L_cache_l1_(1),
        trial_compl_s_U_cache_l1_(1),
        curr_relaxed_compl_x_L_cached_l1_(1),
        curr_relaxed_compl_x_U_cached_l1_(1),
        curr_relaxed_compl_s_L_cached_l1_(1),
        curr_relaxed_compl_s_U_cached_l1_(1),
        curr_complementarity_cache_l1_(6),
        trial_complementarity_cache_l1_(6),
        curr_sigma_x_cache_l1_(1),
        curr_sigma_s_cache_l1_(1),
        dampind_x_L_l1_(NULL),
        dampind_x_U_l1_(NULL),
        dampind_s_L_l1_(NULL),
        dampind_s_U_l1_(NULL),

{
    DBG_START_METH("L1ExactPenaltyRestoCQ::L1ExactPenaltyRestoCQ",
                   dbg_verbosity);
    DBG_ASSERT(IsValid(ip_nlp_) && IsValid(ip_data_));

}
L1ExactPenaltyRestoCQ::~L1ExactPenaltyRestoCQ() noexcept
{ }

void L1ExactPenaltyRestoCQ::RegisterOptions(
        SmartPtr<RegisteredOptions> roptions
        )
{
    roptions->SetRegisteringCategory("Unregistered");
    roptions->AddStringOption4(
            "l1exactpenalty_rho_type",
            "Type of update for the penalty parameter",
            "linear_model",
            "quadratic_model", "check the quadratic model",
            "quadratic_model_no_sigma", "quadratic model without the barrier",
            "linear_model", "use a linear model for predicted reduction",
            "fixed", "fixed rho",
            "Type of update for the penalty."
            );
}

bool L1ExactPenaltyRestoCQ::Initialize1(const Journalist &jnlst,
                                       const OptionsList &options,
                                       const std::string &prefix)
{
    Index rho_int;
    //options.GetEnumValue("l1exactpenalty_rho_type", rho_int, prefix);
//    l1exactpenalty_rho_type_ = RhoUpdateKind(rho_int);
    bool retval = false;

    retval = IpoptCalculatedQuantities::Initialize(jnlst, options, prefix);
    return retval;

}

Number L1ExactPenaltyRestoCQ::curr_f()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_f()",
                   dbg_verbosity)
    Number result;
    SmartPtr<const Vector> x = ip_data_l1_->curr()->x();

    DBG_PRINT_VECTOR(2, "curr_x", *x);
    DBG_PRINT((1, "curr_x tag = %u\n", x->GetTag()));

    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);

    L1ExactPenaltyRestoData& l1data = L1EPRestoData();
    Number rho = L1EPRestoData().CurrentRho();
    sdeps[0] = rho;

    if(!curr_f_cache_l1_.GetCachedResult(result, tdeps, sdeps)){
        if(!trial_f_cache_l1_.GetCachedResult(result, tdeps, sdeps)){
            DBG_PRINT((2, "evaluate current_f\n"));
            result = ip_nlp_l1_->f(*x, rho);
            curr_f_cache_l1_.AddCachedResult(result, tdeps, sdeps);
        }
    }
    DBG_PRINT((1, "result (curr_f) = %e\n", result));
    return result;
}

Number L1ExactPenaltyRestoCQ::trial_f()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::trial_f()",
                   dbg_verbosity);
    Number result;
    SmartPtr<const Vector> x = ip_data_l1_->trial()->x();

    L1ExactPenaltyRestoData& l1data = L1EPRestoData();
    Number rho = l1data.CurrentRho();

    DBG_PRINT_VECTOR(2, "trial_x", *x);
    DBG_PRINT((1, "trial_x tag = %u\n", x-GetTag()));

    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);
    sdeps[0] = rho;

    if( !trial_f_cache_l1_.GetCachedResult(result, tdeps, sdeps)){
        if( !curr_f_cache_l1_.GetCachedResult(result, tdeps, sdeps)){
            DBG_PRINT((2, "evaluate trial_f\n"));
            result = ip_nlp_l1_->f(*x, rho);
            trial_f_cache_l1_.AddCachedResult(result, tdeps, sdeps);
        }
    }
    DBG_PRINT((1, "result (trial_f) = %e\n", result));
    return result;
}

SmartPtr<const Vector> L1ExactPenaltyRestoCQ::curr_grad_f()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_grad_f()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_l1_->curr()->x();

    L1ExactPenaltyRestoData& l1data = L1EPRestoData();
    Number rho = l1data.CurrentRho();

    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);
    sdeps[0] = rho;
    if( !curr_grad_f_cache_l1_.GetCachedResult(result, tdeps, sdeps))
    {
        if( !trial_grad_f_cache_l1_.GetCachedResult(result, tdeps, sdeps))
        {
            result = ip_nlp_l1_->grad_f(*x, rho);
        }
        curr_grad_f_cache_l1_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
}

SmartPtr<const Vector> L1ExactPenaltyRestoCQ::trial_grad_f()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::tial_grad_f",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_l1_->trial()->x();

    L1ExactPenaltyRestoData& l1data = L1EPRestoData();
    Number rho = l1data.CurrentRho();

    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);

    if( !trial_grad_f_cache_l1_.GetCachedResult(result, tdeps, sdeps))
    {
        if( !curr_grad_f_cache_l1_.GetCachedResult(result, tdeps, sdeps))
        {
            result = ip_nlp_l1_->grad_f(*x, rho);
        }
    }
    return result;
}

Number L1ExactPenaltyRestoCQ::CalcBarrierTermL1(
    Number mu,
    const Vector &slack_x_L,
    const Vector &slack_x_U,
    const Vector &slack_s_L,
    const Vector &slack_s_U)
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::CalcBarrierTermL1",
                   dbg_verbosity);
    DBG_START_METH("IpoptCalculatedQuantities::CalcBarrierTermL1",
                   dbg_verbosity);
    //DBG_ASSERT(initialize_called_);

    DBG_PRINT_VECTOR(2, "slack_x_L", slack_x_L);
    DBG_PRINT_VECTOR(2, "slack_x_U", slack_x_U);
    DBG_PRINT_VECTOR(2, "slack_s_L", slack_s_L);
    DBG_PRINT_VECTOR(2, "slack_s_U", slack_s_U);

    Number scale_fact = 1.;
    if( L1EPRestoNlp()->l1exactpenalty_inv_objective_type() )
    {
        scale_fact = 1./L1EPRestoData().CurrentRho();
    }

    Number retval = 0.;
    Number retval2 = 0.;

    const CompoundVector* comp_slack_x_L = static_cast<const CompoundVector*>(&slack_x_L);
    SmartPtr<const Vector> sl_x_only_x_L = comp_slack_x_L->GetComp(0);

    const CompoundVector* comp_slack_x_U = static_cast<const CompoundVector*>(&slack_x_U);
    SmartPtr<const Vector> sl_x_only_x_U = comp_slack_x_U->GetComp(0);

    retval += slack_x_L.SumLogs() - sl_x_only_x_L->SumLogs();
    DBG_PRINT((1, "BarrierTerm after x_L (p + n) = %25.16e\n", retval));
    retval += slack_x_U.SumLogs() - sl_x_only_x_U->SumLogs();
    DBG_PRINT((1, "BarrierTerm after x_U (p + n) = %25.16e\n", retval));
    retval *= -mu;

    retval2 += sl_x_only_x_L->SumLogs();
    DBG_PRINT((1, "BarrierTerm after x_L (x only) = %25.16e\n", retval2));
    retval2 += sl_x_only_x_U->SumLogs();
    DBG_PRINT((1, "BarrierTerm after x_U (x only) = %25.16e\n", retval2));

    retval2 += slack_s_L.SumLogs();
    DBG_PRINT((1, "BarrierTerm after s_L (x only) = %25.16e\n", retval2));
    retval2 += slack_s_U.SumLogs();
    DBG_PRINT((1, "BarrierTerm after s_U (x only) = %25.16e\n", retval2));
    retval2 *= -mu * scale_fact;

    retval += retval2;
    DBG_PRINT((1, "BarrierTerm without damping = %25.16e\n", retval));

    if (kappa_d_l1_ > 0 )
    {
        SmartPtr<const Vector> dampind_x_L;
        SmartPtr<const Vector> dampind_x_U;
        SmartPtr<const Vector> dampind_s_L;
        SmartPtr<const Vector> dampind_s_U;
        ComputeDampingIndicatorsL1(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

        Tmp_x_L_l1().Copy(slack_x_L);
        Tmp_x_L_l1().ElementWiseMultiply(*dampind_x_L);
        retval += kappa_d_l1_ * mu * Tmp_x_L_l1().Asum();
        Tmp_x_U_l1().Copy(slack_x_U);
        Tmp_x_U_l1().ElementWiseMultiply(*dampind_x_U);
        retval += kappa_d_l1_ * mu * Tmp_x_U_l1().Asum();
        Tmp_s_L_l1().Copy(slack_s_L);
        Tmp_s_L_l1().ElementWiseMultiply(*dampind_s_L);
        retval += kappa_d_l1_ * mu * Tmp_s_L_l1().Asum();
        Tmp_s_U_l1().Copy(slack_s_U);
        Tmp_s_U_l1().ElementWiseMultiply(*dampind_s_U);
        retval += kappa_d_l1_ * mu * Tmp_s_U_l1().Asum();
    }
    DBG_PRINT((1, "BarrierTerm with damping = %25.16e\n", retval));

    DBG_ASSERT(IsFiniteNumber(retval));
    return retval;
}

Number L1ExactPenaltyRestoCQ::curr_barrier_obj()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_barrier_obj()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_l1_->curr()->x();
    SmartPtr<const Vector> s = ip_data_l1_->curr()->s();
    DBG_PRINT_VECTOR(2, "curr_x", *x);
    DBG_PRINT_VECTOR(2, "curr_s", *s);
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);

    Number mu = ip_data_l1_->curr_mu();
    Number rho = L1EPRestoData().CurrentRho();
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = rho;

    if( !curr_barrier_obj_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
    {
        if( !trial_barrier_obj_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
        {
            result = curr_f();
            DBG_PRINT((1, "curr_F=%e\n", result));
            result += CalcBarrierTermL1(mu, *curr_slack_x_L(), *curr_slack_x_U(), *curr_slack_s_L(), *curr_slack_s_U());
        }
        curr_barrier_obj_cache_l1_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_ASSERT(IsFiniteNumber(result));

    return result;
}

Number L1ExactPenaltyRestoCQ::trial_barrier_obj()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_barrier_obj()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_l1_->trial()->x();
    SmartPtr<const Vector> s = ip_data_l1_->trial()->s();
    DBG_PRINT_VECTOR(2, "curr_x", *x);
    DBG_PRINT_VECTOR(2, "curr_s", *s);
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);

    Number mu = ip_data_l1_->curr_mu();
    Number rho = L1EPRestoData().CurrentRho();
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = rho;

    if( !trial_barrier_obj_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
    {
        if( !curr_barrier_obj_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
        {
            result = curr_f();
            DBG_PRINT((1, "curr_F=%e\n", result));
            result += CalcBarrierTermL1(mu, *curr_slack_x_L(), *curr_slack_x_U(), *curr_slack_s_L(), *curr_slack_s_U());
        }
        trial_barrier_obj_cache_l1_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_ASSERT(IsFiniteNumber(result));
    return result;
}

SmartPtr<const Vector> L1ExactPenaltyRestoCQ::curr_grad_barrier_obj_x()

{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_grad_barrier_obj_x",
                   dbg_verbosity);
    //DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_l1_->curr()->x();
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    Number mu = ip_data_l1_->curr_mu();
    Number rho = L1EPRestoData().CurrentRho();
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = rho;

    if( !curr_grad_barrier_obj_x_cache_l1_.GetCachedResult(result, tdeps, sdeps))
    {
        SmartPtr<Vector> tmp1 = x->MakeNew();
        tmp1->Copy(*curr_grad_f());

        SmartPtr<Vector> tmp2 = x->MakeNew();
        tmp2->Set(0.0);

        Tmp_x_L_l1().Set(1.0);
        ip_nlp_l1_->Px_L()->AddMSinvZ(-mu, *curr_slack_x_L(), Tmp_x_L_l1(), *tmp2);

        Tmp_x_U_l1().Set(1.0);
        ip_nlp_l1_->Px_U()->AddMSinvZ(mu, *curr_slack_x_U(), Tmp_x_U_l1(), *tmp2);

        if( L1EPRestoNlp()->l1exactpenalty_inv_objective_type())
        {
            CompoundVector* Ctmp2 = static_cast<CompoundVector*>(GetRawPtr(tmp2));
            SmartPtr<Vector> Ctmp2_xonly = Ctmp2->GetCompNonConst(0);
            Ctmp2_xonly->Scal(1/rho);
        }
        tmp1->AddOneVector(1, *tmp2, 1);
        DBG_PRINT_VECTOR(2, "Barrier_Grad_x without damping", *tmp1);

        if( kappa_d_l1_ > 0. )
        {
            SmartPtr<const Vector> dampind_x_L;
            SmartPtr<const Vector> dampind_x_U;
            SmartPtr<const Vector> dampind_s_L;
            SmartPtr<const Vector> dampind_s_U;
            ComputeDampingIndicatorsL1(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

            DBG_PRINT((1, "kappa_d*mu = %e\n", kappa_d_l1_ * mu));
            DBG_PRINT_VECTOR(2, "dampind_x_L", *dampind_x_L);
            ip_nlp_l1_->Px_L()->MultVector(kappa_d_l1_ * mu, *dampind_x_L, 1., *tmp1);
            ip_nlp_l1_->Px_U()->MultVector(-kappa_d_l1_ * mu, *dampind_x_U, 1., *tmp1);
        }

        DBG_PRINT_VECTOR(2, "Barrier_Grad_x with damping", *tmp1);

        result = ConstPtr(tmp1);

        curr_grad_barrier_obj_x_cache_l1_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;

    }

SmartPtr<const Vector> L1ExactPenaltyRestoCQ::curr_grad_barrier_obj_s()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_grad_barrier_obj_s()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> s = ip_data_l1_->curr()->s();
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(s);
    Number mu = ip_data_l1_->curr_mu();
    Number rho = L1EPRestoData().CurrentRho();
    std::vector<Number> sdeps(2);

    sdeps[0] = mu;
    sdeps[1] = rho;

    DBG_PRINT((1, "curr_mu=%e\n", mu));
    DBG_PRINT((1, "curr_rho=%e\n", rho));

    if( !curr_grad_barrier_obj_s_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
    {
        SmartPtr<Vector> tmp1 = s->MakeNew();

        Tmp_s_L_l1().Set(-mu);
        Tmp_s_L_l1().ElementWiseDivide(*curr_slack_s_L());
        ip_nlp_l1_->Pd_L()->MultVector(1., Tmp_s_L_l1(), 0., *tmp1);

        Tmp_s_U_l1().Set(1.);
        ip_nlp_l1_->Pd_U()->AddMSinvZ(mu, *curr_slack_s_U(), Tmp_s_U_l1(), *tmp1);

        // rho part.
        if( L1EPRestoNlp()->l1exactpenalty_inv_objective_type())
        {
            tmp1->Scal(1/rho);
        }

        DBG_PRINT_VECTOR(2, "Barrier_Grad_s without damping", *tmp1);

        // Take care of linear damping terms
        if( kappa_d_l1_ > 0. )
        {
            SmartPtr<const Vector> dampind_x_L;
            SmartPtr<const Vector> dampind_x_U;
            SmartPtr<const Vector> dampind_s_L;
            SmartPtr<const Vector> dampind_s_U;
            ComputeDampingIndicatorsL1(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

            DBG_PRINT((1, "kappa_d*mu = %e\n", kappa_d_l1_ * mu));
            DBG_PRINT_VECTOR(2, "dampind_s_L", *dampind_s_L);
            DBG_PRINT_VECTOR(2, "dampind_s_U", *dampind_s_U);
            ip_nlp_l1_->Pd_L()->MultVector(kappa_d_l1_ * mu, *dampind_s_L, 1., *tmp1);
            ip_nlp_l1_->Pd_U()->MultVector(-kappa_d_l1_ * mu, *dampind_s_U, 1., *tmp1);
        }

        DBG_PRINT_VECTOR(2, "Barrier_Grad_s with damping", *tmp1);

        result = ConstPtr(tmp1);

        curr_grad_barrier_obj_s_cache_l1_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
}

SmartPtr<const Vector> L1ExactPenaltyRestoCQ::grad_kappa_times_damping_x()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::grad_kappa_times_damping_x()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_l1_->curr()->x();

    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(ip_nlp_l1_->Px_L());
    tdeps[1] = GetRawPtr(ip_nlp_l1_->Px_U());
    std::vector<Number> sdeps(1);
    sdeps[0] = L1EPRestoData().CurrentRho();
    if( !grad_kappa_times_damping_x_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
    {
        SmartPtr<Vector> tmp1 = x->MakeNew();
        if( kappa_d_l1_ > 0. )
        {
            SmartPtr<const Vector> dampind_x_L;
            SmartPtr<const Vector> dampind_x_U;
            SmartPtr<const Vector> dampind_s_L;
            SmartPtr<const Vector> dampind_s_U;
            ComputeDampingIndicatorsL1(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

            ip_nlp_l1_->Px_L()->MultVector(kappa_d_l1_, *dampind_x_L, 0., *tmp1);
            ip_nlp_l1_->Px_U()->MultVector(-kappa_d_l1_, *dampind_x_U, 1., *tmp1);
        }
        else
        {
            tmp1->Set(0.);
        }
        result = ConstPtr(tmp1);

        grad_kappa_times_damping_x_cache_l1_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
}

SmartPtr<const Vector> L1ExactPenaltyRestoCQ::grad_kappa_times_damping_s()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::grad_kappa_times_damping_s()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> s = ip_data_l1_->curr()->s();

    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(ip_nlp_l1_->Pd_L());
    tdeps[1] = GetRawPtr(ip_nlp_l1_->Pd_U());
    std::vector<Number> sdeps(1);
    sdeps[0] = L1EPRestoData().CurrentRho();
    if( !grad_kappa_times_damping_s_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
    {
        SmartPtr<Vector> tmp1 = s->MakeNew();
        if( kappa_d_l1_ > 0. )
        {
            SmartPtr<const Vector> dampind_x_L;
            SmartPtr<const Vector> dampind_x_U;
            SmartPtr<const Vector> dampind_s_L;
            SmartPtr<const Vector> dampind_s_U;
            ComputeDampingIndicatorsL1(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

            ip_nlp_l1_->Pd_L()->MultVector(kappa_d_l1_, *dampind_s_L, 0., *tmp1);
            ip_nlp_l1_->Pd_U()->MultVector(-kappa_d_l1_, *dampind_s_U, 1., *tmp1);
        }
        else
        {
            tmp1->Set(0.);
        }
        result = ConstPtr(tmp1);

        grad_kappa_times_damping_s_cache_l1_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
}

void L1ExactPenaltyRestoCQ::ComputeDampingIndicatorsL1(
        SmartPtr<const Vector> &dampind_x_L,
        SmartPtr<const Vector> &dampind_x_U,
        SmartPtr<const Vector> &dampind_s_L,
        SmartPtr<const Vector> &dampind_s_U)
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::ComputeDampingIndicatorsL1",
                   dbg_verbosity);

    Number scale = 1.0;
    SmartPtr<const Vector> result;

    if( L1EPRestoNlp()->l1exactpenalty_inv_objective_type())
    {
        scale = 1./L1ExactPenaltyRestoData().CurrentRho();
        DBG_PRINT((2, "Inverse scaling! \n"));
    }
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = NULL ;

    std::vector<Number> sdeps(1);
    sdeps[0] = L1ExactPenaltyRestoData().CurrentRho();
    dampind_s_U_cache_l1_.GetCachedResult(result, tdeps, sdeps);

    DBG_PRINT((2, "scale_damping_l1_ = %e\n", scale));
    // What if rho changes?: We add caches.
    //if( IsNull(dampind_x_L_l1_))
    if(!dampind_x_L_cache_l1_.GetCachedResult(dampind_x_L, tdeps, sdeps)) {
        Tmp_x_L_l1().Set(1.0);
        ip_nlp_l1_->Px_L()->MultVector(1.0, Tmp_x_L_l1(), 0.0, Tmp_x_l1());
        Tmp_x_U_l1().Set(1.0);
        ip_nlp_l1_->Px_U()->MultVector(-1.0, Tmp_x_U_l1(), 1.0, Tmp_x_l1());

        dampind_x_L_l1_ = ip_nlp_l1_->x_L()->MakeNew();
        ip_nlp_l1_->Px_L()->TransMultVector(1.0, Tmp_x_l1(), 0.0,
                                            *dampind_x_L_l1_);

        dampind_x_U_l1_ = ip_nlp_l1_->x_U()->MakeNew();
        ip_nlp_l1_->Px_U()->TransMultVector(-1.0, Tmp_x_l1(), 0.0,
                                            *dampind_x_U_l1_);

        // Now for s
        Tmp_s_L_l1().Set(1.0);
        ip_nlp_l1_->Pd_L()->MultVector(1.0, Tmp_s_L_l1(), 0.0, Tmp_s_l1());
        Tmp_s_U_l1().Set(1.0);
        ip_nlp_l1_->Pd_U()->MultVector(-1.0, Tmp_s_U_l1(), 1.0, Tmp_s_l1());

        dampind_s_L_l1_ = ip_nlp_l1_->d_L()->MakeNew();
        ip_nlp_l1_->Pd_L()->TransMultVector(1.0, Tmp_s_l1(), 0.0,
                                            *dampind_s_L_l1_);

        dampind_s_U_l1_ = ip_nlp_l1_->d_U()->MakeNew();
        ip_nlp_l1_->Pd_U()->TransMultVector(-1.0, Tmp_s_l1(), 0.0,
                                            *dampind_s_U_l1_);

        DBG_PRINT_VECTOR(2, "dampind_x_L_l1_", *dampind_x_L_l1_);
        DBG_PRINT_VECTOR(2, "dampind_x_U_l1_", *dampind_x_U_l1_);
        DBG_PRINT_VECTOR(2, "dampind_s_L_l1_", *dampind_s_L_l1_);
        DBG_PRINT_VECTOR(2, "dampind_s_U_l1_", *dampind_s_U_l1_);

        // Scale the relevant parts if necessary.
        CompoundVector *C_dampind_x_L_ = static_cast<CompoundVector *>(GetRawPtr(
                dampind_x_L_l1_));
        SmartPtr <Vector> tmp_xonly_x_L_ = C_dampind_x_L_->GetCompNonConst(0);

        CompoundVector *C_dampind_x_U_ = static_cast<CompoundVector *>(GetRawPtr(
                dampind_x_U_l1_));
        SmartPtr <Vector> tmp_xonly_x_U_ = C_dampind_x_U_->GetCompNonConst(0);

        // We got to introduce the 1/rho here for the damping terms.
        // For both xL and xU slacks and sL, and sU.
        // The latter are not compound vectors.

        tmp_xonly_x_L_->Scal(scale);
        tmp_xonly_x_U_->Scal(scale);
        dampind_s_L_l1_->Scal(scale);
        dampind_s_U_l1_->Scal(scale);

        DBG_PRINT_VECTOR(2, "(scaled)dampind_x_L_l1_", *dampind_x_L_l1_);
        DBG_PRINT_VECTOR(2, "(scaled)dampind_x_U_l1_", *dampind_x_U_l1_);
        DBG_PRINT_VECTOR(2, "(scaled)dampind_s_L_l1_", *dampind_s_L_l1_);
        DBG_PRINT_VECTOR(2, "(scaled)dampind_s_U_l1_", *dampind_s_U_l1_);

        dampind_x_L = ConstPtr(dampind_x_L_l1_);
        dampind_x_U = ConstPtr(dampind_x_U_l1_);
        dampind_s_L = ConstPtr(dampind_s_L_l1_);
        dampind_s_U = ConstPtr(dampind_s_U_l1_);

        dampind_x_L_cache_l1_.AddCachedResult(dampind_x_L, tdeps, sdeps);
        dampind_x_U_cache_l1_.AddCachedResult(dampind_x_U, tdeps, sdeps);
        dampind_s_L_cache_l1_.AddCachedResult(dampind_s_L, tdeps, sdeps);
        dampind_s_U_cache_l1_.AddCachedResult(dampind_s_U, tdeps, sdeps);

    } else {
        dampind_x_L_cache_l1_.GetCachedResult(dampind_x_L, tdeps, sdeps);
        dampind_x_U_cache_l1_.GetCachedResult(dampind_x_U, tdeps, sdeps);
        dampind_s_L_cache_l1_.GetCachedResult(dampind_s_L, tdeps, sdeps);
        dampind_s_U_cache_l1_.GetCachedResult(dampind_s_U, tdeps, sdeps);
    }


}

SmartPtr<const SymMatrix> L1ExactPenaltyRestoCQ::curr_exact_hessian()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_exact_hessian()",
                   dbg_verbosity);

    SmartPtr<const SymMatrix> result;
    SmartPtr<const Vector> x = ip_data_l1_->curr()->x();
    SmartPtr<const Vector> y_c = ip_data_l1_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_l1_->curr()->y_d();
    std::vector<const TaggedObject*> tdeps(3);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(y_c);
    tdeps[2] = GetRawPtr(y_d);

    Number rho = L1EPRestoData().CurrentRho();
    std::vector<Number> sdeps(1);
    sdeps[0] = rho;

    if( !curr_exact_hessian_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
    {
        result = ip_nlp_l1_->h(*x, 1.0, *y_c, *y_d, rho);
        curr_exact_hessian_cache_l1_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
}

SmartPtr<const Vector> L1ExactPenaltyRestoCQ::curr_grad_lag_x()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_grad_lag_x",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_l1_->curr()->x();
    SmartPtr<const Vector> y_c = ip_data_l1_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_l1_->curr()->y_d();
    SmartPtr<const Vector> z_L = ip_data_l1_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_l1_->curr()->z_U();

    Number rho = L1EPRestoData().CurrentRho();

    std::vector<const TaggedObject*> tdeps(5);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(y_c);
    tdeps[2] = GetRawPtr(y_d);
    tdeps[3] = GetRawPtr(z_L);
    tdeps[4] = GetRawPtr(z_U);

    std::vector<Number> sdeps(1);
    sdeps[0] = rho;

    if( !curr_grad_lag_x_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
    {
        if( !trial_grad_lag_x_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
        {
            SmartPtr<Vector> tmp = x->MakeNew();
            DBG_PRINT_VECTOR(2, "curr_grad_f", *curr_grad_f());
            tmp->Copy(*curr_grad_f());
            tmp->AddTwoVectors(1., *curr_jac_cT_times_curr_y_c(), 1., *curr_jac_dT_times_curr_y_d(), 1.);
            DBG_PRINT_VECTOR(2, "jac_cT*y_c", *curr_jac_cT_times_curr_y_c());
            DBG_PRINT_VECTOR(2, "jac_dT*y_d", *curr_jac_dT_times_curr_y_d());

            if( L1EPRestoNlp()->l1exactpenalty_inv_objective_type())
            {   // Take care of the (-z_L + z_U)/rho part.
                SmartPtr<Vector> tmpzL = z_L->MakeNewCopy();
                SmartPtr<Vector> tmpzU = z_U->MakeNewCopy();
                CompoundVector* C_z_L = static_cast<CompoundVector*>(GetRawPtr(tmpzL));
                DBG_ASSERT(dynamic_cast<CompundVector*>(GetRawPtr(tmpzL)));
                SmartPtr<Vector> z_L_x_only = C_z_L->GetCompNonConst(0);
                C_z_L->Scal(1/rho);
                CompoundVector* C_z_U = static_cast<CompoundVector*>(GetRawPtr(tmpzU));
                DBG_ASSERT(dynamic_cast<CompundVector*>(GetRawPtr(tmpzU)));
                SmartPtr<Vector> z_U_x_only = C_z_U->GetCompNonConst(0);
                C_z_U->Scal(1/rho);
                ip_nlp_l1_->Px_L()->MultVector(-1., *tmpzL, 1., *tmp);
                ip_nlp_l1_->Px_U()->MultVector(1., *tmpzU, 1., *tmp);
            }
            else {
                ip_nlp_l1_->Px_L()->MultVector(-1., *z_L, 1., *tmp);
                ip_nlp_l1_->Px_U()->MultVector(1., *z_U, 1., *tmp);
            }

            result = ConstPtr(tmp);
        }
        curr_grad_lag_x_cache_l1_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;

}

SmartPtr<const Vector> L1ExactPenaltyRestoCQ::trial_grad_lag_x()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::trial_grad_lag_x",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_l1_->trial()->x();
    SmartPtr<const Vector> y_c = ip_data_l1_->trial()->y_c();
    SmartPtr<const Vector> y_d = ip_data_l1_->trial()->y_d();
    SmartPtr<const Vector> z_L = ip_data_l1_->trial()->z_L();
    SmartPtr<const Vector> z_U = ip_data_l1_->trial()->z_U();

    Number rho = L1EPRestoData().CurrentRho();

    std::vector<const TaggedObject*> tdeps(5);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(y_c);
    tdeps[2] = GetRawPtr(y_d);
    tdeps[3] = GetRawPtr(z_L);
    tdeps[4] = GetRawPtr(z_U);

    std::vector<Number> sdeps(1);
    sdeps[0] = rho;

    if( !trial_grad_lag_x_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
    {
        if( !curr_grad_lag_x_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
        {
            SmartPtr<Vector> tmp = x->MakeNew();
            DBG_PRINT_VECTOR(2, "curr_grad_f", *curr_grad_f());
            tmp->Copy(*curr_grad_f());
            tmp->AddTwoVectors(1., *curr_jac_cT_times_curr_y_c(), 1., *curr_jac_dT_times_curr_y_d(), 1.);
            DBG_PRINT_VECTOR(2, "jac_cT*y_c", *curr_jac_cT_times_curr_y_c());
            DBG_PRINT_VECTOR(2, "jac_dT*y_d", *curr_jac_dT_times_curr_y_d());

            if( L1EPRestoNlp()->l1exactpenalty_inv_objective_type())
            {   // Take care of the (-z_L + z_U)/rho part.
                SmartPtr<Vector> tmpzL = z_L->MakeNewCopy();
                SmartPtr<Vector> tmpzU = z_U->MakeNewCopy();
                CompoundVector* C_z_L = static_cast<CompoundVector*>(GetRawPtr(tmpzL));
                DBG_ASSERT(dynamic_cast<CompundVector*>(GetRawPtr(tmpzL)));
                SmartPtr<Vector> z_L_x_only = C_z_L->GetCompNonConst(0);
                C_z_L->Scal(1/rho);
                CompoundVector* C_z_U = static_cast<CompoundVector*>(GetRawPtr(tmpzU));
                DBG_ASSERT(dynamic_cast<CompundVector*>(GetRawPtr(tmpzU)));
                SmartPtr<Vector> z_U_x_only = C_z_U->GetCompNonConst(0);
                C_z_U->Scal(1/rho);
                ip_nlp_l1_->Px_L()->MultVector(-1., *tmpzL, 1., *tmp);
                ip_nlp_l1_->Px_U()->MultVector(1., *tmpzU, 1., *tmp);
            }
            else {
                ip_nlp_l1_->Px_L()->MultVector(-1., *z_L, 1., *tmp);
                ip_nlp_l1_->Px_U()->MultVector(1., *z_U, 1., *tmp);
            }

            result = ConstPtr(tmp);
        }
        trial_grad_lag_x_cache_l1_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;

}

SmartPtr<const Vector> L1ExactPenaltyRestoCQ::curr_grad_lag_s()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_grad_lag_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> y_d = ip_data_l1_->curr()->y_d();
    SmartPtr<const Vector> v_L = ip_data_l1_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_l1_->curr()->v_U();

    std::vector<const TaggedObject*> deps(3);
    deps[0] = GetRawPtr(y_d);
    deps[1] = GetRawPtr(v_L);
    deps[2] = GetRawPtr(v_U);

    Number rho = L1EPRestoData().CurrentRho();
    std::vector<Number> sdeps(1);
    sdeps[0] = rho;

    if( !curr_grad_lag_s_cache_l1_.GetCachedResult(result, deps, sdeps) )
    {
        if( !trial_grad_lag_s_cache_l1_.GetCachedResult(result, deps, sdeps) )
        {
            SmartPtr<Vector> tmp = y_d->MakeNew();
            ip_nlp_l1_->Pd_U()->MultVector(1., *v_U, 0., *tmp);
            ip_nlp_l1_->Pd_L()->MultVector(-1., *v_L, 1., *tmp);
            if(L1EPRestoNlp()->l1exactpenalty_inv_objective_type())
            {
                tmp->Scal(1/rho);
            }
            tmp->Axpy(-1., *y_d);
            result = ConstPtr(tmp);
        }
        curr_grad_lag_s_cache_l1_.AddCachedResult(result, deps, sdeps);
    }

    return result;
}

SmartPtr<const Vector> L1ExactPenaltyRestoCQ::trial_grad_lag_s()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::trial_grad_lag_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> y_d = ip_data_l1_->trial()->y_d();
    SmartPtr<const Vector> v_L = ip_data_l1_->trial()->v_L();
    SmartPtr<const Vector> v_U = ip_data_l1_->trial()->v_U();

    std::vector<const TaggedObject*> deps(3);
    deps[0] = GetRawPtr(y_d);
    deps[1] = GetRawPtr(v_L);
    deps[2] = GetRawPtr(v_U);

    Number rho = L1EPRestoData().CurrentRho();
    std::vector<Number> sdeps(1);
    sdeps[0] = rho;

    if( !trial_grad_lag_s_cache_l1_.GetCachedResult(result, deps, sdeps) )
    {
        if( !curr_grad_lag_s_cache_l1_.GetCachedResult(result, deps, sdeps) )
        {
            SmartPtr<Vector> tmp = y_d->MakeNew();
            ip_nlp_l1_->Pd_U()->MultVector(1., *v_U, 0., *tmp);
            ip_nlp_l1_->Pd_L()->MultVector(-1., *v_L, 1., *tmp);
            if(L1EPRestoNlp()->l1exactpenalty_inv_objective_type())
            {
                tmp->Scal(1/rho);
            }
            tmp->Axpy(-1., *y_d);
            result = ConstPtr(tmp);
        }
        trial_grad_lag_s_cache_l1_.AddCachedResult(result, deps, sdeps);
    }

    return result;
}


SmartPtr<const Vector> L1ExactPenaltyRestoCQ::curr_grad_lag_with_damping_x()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_grad_lag_with_damping_x()",
                   dbg_verbosity);

    /* If no damping is used, just return the gradient of the regular
     Lagrangian function */
    if( kappa_d_l1_ == 0. )
    {
        return curr_grad_lag_x();
    }

    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_l1_->curr()->x();
    SmartPtr<const Vector> y_c = ip_data_l1_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_l1_->curr()->y_d();
    SmartPtr<const Vector> z_L = ip_data_l1_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_l1_->curr()->z_U();
    Number mu = ip_data_l1_->curr_mu();
    Number rho = L1EPRestoData().CurrentRho();
    std::vector<const TaggedObject*> deps(5);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(y_c);
    deps[2] = GetRawPtr(y_d);
    deps[3] = GetRawPtr(z_L);
    deps[4] = GetRawPtr(z_U);
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = rho;

    if( !curr_grad_lag_with_damping_x_cache_l1_.GetCachedResult(result, deps, sdeps) )
    {
        SmartPtr<Vector> tmp = x->MakeNew();
        tmp->Copy(*curr_grad_lag_x());

        SmartPtr<const Vector> dampind_x_L;
        SmartPtr<const Vector> dampind_x_U;
        SmartPtr<const Vector> dampind_s_L;
        SmartPtr<const Vector> dampind_s_U;
        ComputeDampingIndicatorsL1(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

        ip_nlp_l1_->Px_L()->MultVector(kappa_d_l1_ * mu, *dampind_x_L, 1., *tmp);
        ip_nlp_l1_->Px_U()->MultVector(-kappa_d_l1_ * mu, *dampind_x_U, 1., *tmp);

        result = ConstPtr(tmp);
        curr_grad_lag_with_damping_x_cache_l1_.AddCachedResult(result, deps, sdeps);
    }

    return result;
}


SmartPtr<const Vector> L1ExactPenaltyRestoCQ::curr_grad_lag_with_damping_s()
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::curr_grad_lag_with_damping_s()",
                   dbg_verbosity);

    /* If no damping is used, just return the gradient of the regular
     Lagrangian function */
    if( kappa_d_l1_ == 0. )
    {
        return curr_grad_lag_s();
    }

    SmartPtr<const Vector> result;

    SmartPtr<const Vector> y_d = ip_data_l1_->curr()->y_d();
    SmartPtr<const Vector> v_L = ip_data_l1_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_l1_->curr()->v_U();
    std::vector<const TaggedObject*> deps(3);
    deps[0] = GetRawPtr(y_d);
    deps[1] = GetRawPtr(v_L);
    deps[2] = GetRawPtr(v_U);

    Number mu = ip_data_l1_->curr_mu();
    Number rho = L1EPRestoData().CurrentRho();
    std::vector<Number> sdeps(2);
    sdeps[0] = mu;
    sdeps[1] = rho;

    if( !curr_grad_lag_with_damping_s_cache_l1_.GetCachedResult(result, deps, sdeps) )
    {
        SmartPtr<Vector> tmp = y_d->MakeNew();
        tmp->Copy(*curr_grad_lag_s());

        SmartPtr<const Vector> dampind_x_L;
        SmartPtr<const Vector> dampind_x_U;
        SmartPtr<const Vector> dampind_s_L;
        SmartPtr<const Vector> dampind_s_U;
        ComputeDampingIndicatorsL1(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

        ip_nlp_l1_->Pd_L()->MultVector(kappa_d_l1_ * mu, *dampind_s_L, 1., *tmp);
        ip_nlp_l1_->Pd_U()->MultVector(-kappa_d_l1_ * mu, *dampind_s_U, 1., *tmp);

        result = ConstPtr(tmp);
        curr_grad_lag_with_damping_s_cache_l1_.AddCachedResult(result, deps, sdeps);
    }

    return result;
}



Vector& L1ExactPenaltyRestoCQ::Tmp_x_l1()
{
    if( !IsValid(tmp_x_l1_) )
    {
        tmp_x_l1_ = ip_data_l1_->curr()->x()->MakeNew();
    }
    return *tmp_x_l1_;
}

Vector& L1ExactPenaltyRestoCQ::Tmp_s_l1()
{
    if( !IsValid(tmp_s_l1_) )
    {
        tmp_s_l1_ = ip_data_l1_->curr()->s()->MakeNew();
    }
    return *tmp_s_l1_;
}

Vector& L1ExactPenaltyRestoCQ::Tmp_c_l1()
{
    if( !IsValid(tmp_c_l1_) )
    {
        tmp_c_l1_ = ip_data_l1_->curr()->y_c()->MakeNew();
    }
    return *tmp_c_l1_;
}

Vector& L1ExactPenaltyRestoCQ::Tmp_d_l1()
{
    if( !IsValid(tmp_d_l1_) )
    {
        tmp_d_l1_ = ip_data_l1_->curr()->y_d()->MakeNew();
    }
    return *tmp_d_l1_;
}

Vector& L1ExactPenaltyRestoCQ::Tmp_x_L_l1()
{
    if( !IsValid(tmp_x_L_l1_) )
    {
        tmp_x_L_l1_ = ip_nlp_l1_->x_L()->MakeNew();
    }
    return *tmp_x_L_l1_;
}

Vector& L1ExactPenaltyRestoCQ::Tmp_x_U_l1()
{
    if( !IsValid(tmp_x_U_l1_) )
    {
        tmp_x_U_l1_ = ip_nlp_l1_->x_U()->MakeNew();
    }
    return *tmp_x_U_l1_;
}

Vector& L1ExactPenaltyRestoCQ::Tmp_s_L_l1()
{
    if( !IsValid(tmp_s_L_l1_) )
    {
        tmp_s_L_l1_ = ip_nlp_l1_->d_L()->MakeNew();
    }
    return *tmp_s_L_l1_;
}

Vector& L1ExactPenaltyRestoCQ::Tmp_s_U_l1()
{
    if( !IsValid(tmp_s_U_l1_) )
    {
        tmp_s_U_l1_ = ip_nlp_l1_->d_U()->MakeNew();
    }
    return *tmp_s_U_l1_;
}

}