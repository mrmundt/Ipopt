//
// Created by David on 10/12/2020.
//

#include "IpL1ExactPenaltyRestoCalculatedQuantities.hpp"

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
        dampind_x_L_l1_(nullptr),
        dampind_x_U_l1_(nullptr),
        dampind_s_L_l1_(nullptr),
        dampind_s_U_l1_(nullptr),
        curr_exact_hessian_cache_l1_(1)   ,
        curr_grad_lag_x_cache_l1_(1),
        trial_grad_lag_x_cache_l1_(1),
        curr_grad_lag_s_cache_l1_(1),
        trial_grad_lag_s_cache_l1_(1),
        curr_grad_lag_with_damping_x_cache_l1_(1),
        curr_grad_lag_with_damping_s_cache_l1_(1)
{
    DBG_START_METH("L1ExactPenaltyRestoCQ::L1ExactPenaltyRestoCQ",
                   dbg_verbosity);
    DBG_ASSERT(IsValid(ip_nlp_) && IsValid(ip_data_));
}
L1ExactPenaltyRestoCQ::~L1ExactPenaltyRestoCQ()
= default;

bool L1ExactPenaltyRestoCQ::Initialize(
        const Journalist&  jnlst,
        const OptionsList& options,
        const std::string& prefix
)
{
    Index enum_int;

    //options.GetNumericValue("s_max", s_max_, prefix);
    options.GetNumericValue("kappa_d", kappa_d_l1_, prefix);
    //options.GetNumericValue("slack_move", slack_move_, prefix);
    options.GetEnumValue("constraint_violation_norm_type", enum_int, prefix);
    //constr_viol_normtype_ = ENormType(enum_int);
    // The following option is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure", warm_start_same_structure_l1_, prefix);
    //options.GetNumericValue("mu_target", mu_target_, prefix);
    jnlst_ = &jnlst;
    if( !warm_start_same_structure_l1_ )
    {
        dampind_x_L_l1_ = nullptr;
        dampind_x_U_l1_ = nullptr;
        dampind_s_L_l1_ = nullptr;
        dampind_s_U_l1_ = nullptr;

        tmp_x_l1_ = nullptr;
        tmp_s_l1_ = nullptr;
        //tmp_c_l1_ = nullptr;
        //tmp_d_l1_ = nullptr;
        tmp_x_L_l1_ = nullptr;
        tmp_x_U_l1_ = nullptr;
        tmp_s_L_l1_ = nullptr;
        tmp_s_U_l1_ = nullptr;
    }


    //initialize_called_ = true;

    bool retval;
    //if( IsValid(add_cq_) )
    //{
    //    retval = add_cq_->Initialize(jnlst, options, prefix);
    //}
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
    Number rho = L1EPRestoData().GetCurrentRho();
    sdeps[0] = rho;
    if(!curr_f_cache_l1_.GetCachedResult(result, tdeps, sdeps)){
        if(!trial_f_cache_l1_.GetCachedResult(result, tdeps, sdeps)){
            DBG_PRINT((2, "evaluate current_f\n"));
            result = L1EPRestoNlp()->f(*x, rho);
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
    Number rho = l1data.GetCurrentRho();
    DBG_PRINT_VECTOR(2, "trial_x", *x);
    DBG_PRINT((1, "trial_x tag = %u\n", x-GetTag()));

    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);
    sdeps[0] = rho;

    if( !trial_f_cache_l1_.GetCachedResult(result, tdeps, sdeps)){
        if( !curr_f_cache_l1_.GetCachedResult(result, tdeps, sdeps)){
            DBG_PRINT((2, "evaluate trial_f\n"));
            result = L1EPRestoNlp()->f(*x, rho);
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
    Number rho = l1data.GetCurrentRho();

    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);
    sdeps[0] = rho;
    if( !curr_grad_f_cache_l1_.GetCachedResult(result, tdeps, sdeps))
    {
        if( !trial_grad_f_cache_l1_.GetCachedResult(result, tdeps, sdeps))
        {
            result = L1EPRestoNlp()->grad_f(*x, rho);
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
    Number rho = l1data.GetCurrentRho();

    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);

    if( !trial_grad_f_cache_l1_.GetCachedResult(result, tdeps, sdeps))
    {
        if( !curr_grad_f_cache_l1_.GetCachedResult(result, tdeps, sdeps))
        {
            result = L1EPRestoNlp()->grad_f(*x, rho);
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
    if(L1EPRestoNlp()->l1_epr_inv_objective_type() )
    {
        scale_fact = 1./ L1EPRestoData().GetCurrentRho();
    }
    Number retval = 0.;
    Number retval2 = 0.;

    auto comp_slack_x_L = dynamic_cast<const CompoundVector*>(&slack_x_L);
    SmartPtr<const Vector> sl_x_only_x_L = comp_slack_x_L->GetComp(0);

    auto comp_slack_x_U = dynamic_cast<const CompoundVector*>(&slack_x_U);
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
    Number rho = L1EPRestoData().GetCurrentRho();
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
    Number rho = L1EPRestoData().GetCurrentRho();
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
    Number rho = L1EPRestoData().GetCurrentRho();
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
        L1EPRestoNlp()->Px_L()->AddMSinvZ(-mu, *curr_slack_x_L(), Tmp_x_L_l1(), *tmp2);

        Tmp_x_U_l1().Set(1.0);
        L1EPRestoNlp()->Px_U()->AddMSinvZ(mu, *curr_slack_x_U(), Tmp_x_U_l1(), *tmp2);

        if(L1EPRestoNlp()->l1_epr_inv_objective_type())
        { // only do it once.
            auto Ctmp2 = dynamic_cast<CompoundVector*>(GetRawPtr(tmp2));
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
            L1EPRestoNlp()->Px_L()->MultVector(kappa_d_l1_ * mu, *dampind_x_L, 1., *tmp1);
            L1EPRestoNlp()->Px_U()->MultVector(-kappa_d_l1_ * mu, *dampind_x_U, 1., *tmp1);
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
    //DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> s = ip_data_l1_->curr()->s();
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(s);
    Number mu = ip_data_l1_->curr_mu();
    Number rho = L1EPRestoData().GetCurrentRho();
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
        L1EPRestoNlp()->Pd_L()->MultVector(1., Tmp_s_L_l1(), 0., *tmp1);

        Tmp_s_U_l1().Set(1.);
        L1EPRestoNlp()->Pd_U()->AddMSinvZ(mu, *curr_slack_s_U(), Tmp_s_U_l1(), *tmp1);

        // rho part.
        if(L1EPRestoNlp()->l1_epr_inv_objective_type())
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
            L1EPRestoNlp()->Pd_L()->MultVector(kappa_d_l1_ * mu, *dampind_s_L, 1., *tmp1);
            L1EPRestoNlp()->Pd_U()->MultVector(-kappa_d_l1_ * mu, *dampind_s_U, 1., *tmp1);
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
    //DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_l1_->curr()->x();

    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(L1EPRestoNlp()->Px_L());
    tdeps[1] = GetRawPtr(L1EPRestoNlp()->Px_U());
    std::vector<Number> sdeps(1);
    sdeps[0] = L1EPRestoData().GetCurrentRho();
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

            L1EPRestoNlp()->Px_L()->MultVector(kappa_d_l1_, *dampind_x_L, 0., *tmp1);
            L1EPRestoNlp()->Px_U()->MultVector(-kappa_d_l1_, *dampind_x_U, 1., *tmp1);
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
    //DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> s = ip_data_l1_->curr()->s();

    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(L1EPRestoNlp()->Pd_L());
    tdeps[1] = GetRawPtr(L1EPRestoNlp()->Pd_U());
    std::vector<Number> sdeps(1);
    sdeps[0] = L1EPRestoData().GetCurrentRho();
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

            L1EPRestoNlp()->Pd_L()->MultVector(kappa_d_l1_, *dampind_s_L, 0., *tmp1);
            L1EPRestoNlp()->Pd_U()->MultVector(-kappa_d_l1_, *dampind_s_U, 1., *tmp1);
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

    Number scale = 1.;
    SmartPtr<const Vector> result;
    if(L1EPRestoNlp()->l1_epr_inv_objective_type())
    {

        scale = 1./ L1EPRestoData().GetCurrentRho();
        DBG_PRINT((2, "Inverse scaling! \n"));
    }
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = nullptr ;

    std::vector<Number> sdeps(1);
    sdeps[0] = L1EPRestoData().GetCurrentRho();
    dampind_s_U_cache_l1_.GetCachedResult(result, tdeps, sdeps);

    DBG_PRINT((2, "scale_damping_l1_ = %e\n", scale));
    // What if rho changes?: We add caches.
    //if( IsNull(dampind_x_L_l1_))
    if(!dampind_x_L_cache_l1_.GetCachedResult(dampind_x_L, tdeps, sdeps)) {
        Tmp_x_L_l1().Set(1.0);
        L1EPRestoNlp()->Px_L()->MultVector(1.0, Tmp_x_L_l1(), 0.0, Tmp_x_l1());
        Tmp_x_U_l1().Set(1.0);
        L1EPRestoNlp()->Px_U()->MultVector(-1.0, Tmp_x_U_l1(), 1.0, Tmp_x_l1());

        dampind_x_L_l1_ = L1EPRestoNlp()->x_L()->MakeNew();
        L1EPRestoNlp()->Px_L()->TransMultVector(1.0, Tmp_x_l1(), 0.0,
                                            *dampind_x_L_l1_);

        dampind_x_U_l1_ = L1EPRestoNlp()->x_U()->MakeNew();
        L1EPRestoNlp()->Px_U()->TransMultVector(-1.0, Tmp_x_l1(), 0.0,
                                            *dampind_x_U_l1_);

        // Now for s
        Tmp_s_L_l1().Set(1.0);
        L1EPRestoNlp()->Pd_L()->MultVector(1.0, Tmp_s_L_l1(), 0.0, Tmp_s_l1());
        Tmp_s_U_l1().Set(1.0);
        L1EPRestoNlp()->Pd_U()->MultVector(-1.0, Tmp_s_U_l1(), 1.0, Tmp_s_l1());

        dampind_s_L_l1_ = L1EPRestoNlp()->d_L()->MakeNew();
        L1EPRestoNlp()->Pd_L()->TransMultVector(1.0, Tmp_s_l1(), 0.0,
                                            *dampind_s_L_l1_);

        dampind_s_U_l1_ = L1EPRestoNlp()->d_U()->MakeNew();
        L1EPRestoNlp()->Pd_U()->TransMultVector(-1.0, Tmp_s_l1(), 0.0,
                                            *dampind_s_U_l1_);

        DBG_PRINT_VECTOR(2, "dampind_x_L_l1_", *dampind_x_L_l1_);
        DBG_PRINT_VECTOR(2, "dampind_x_U_l1_", *dampind_x_U_l1_);
        DBG_PRINT_VECTOR(2, "dampind_s_L_l1_", *dampind_s_L_l1_);
        DBG_PRINT_VECTOR(2, "dampind_s_U_l1_", *dampind_s_U_l1_);

        // Scale the relevant parts if necessary.
        auto C_dampind_x_L_ = dynamic_cast<CompoundVector *>(GetRawPtr(
                dampind_x_L_l1_));
        SmartPtr <Vector> tmp_xonly_x_L_ = C_dampind_x_L_->GetCompNonConst(0);

        auto C_dampind_x_U_ = dynamic_cast<CompoundVector *>(GetRawPtr(
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

    Number rho = L1EPRestoData().GetCurrentRho();
    std::vector<Number> sdeps(1);
    sdeps[0] = rho;

    if( !curr_exact_hessian_cache_l1_.GetCachedResult(result, tdeps, sdeps) )
    {
        result = L1EPRestoNlp()->h(*x, 1.0, *y_c, *y_d, rho);
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

    Number rho = L1EPRestoData().GetCurrentRho();
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

            if(L1EPRestoNlp()->l1_epr_inv_objective_type())
            {   // Take care of the (-z_L + z_U)/rho part.
                SmartPtr<Vector> tmpzL = z_L->MakeNewCopy();
                auto C_z_L = dynamic_cast<CompoundVector*>(GetRawPtr(tmpzL));
                DBG_ASSERT(dynamic_cast<CompundVector*>(GetRawPtr(tmpzL)));
                SmartPtr<Vector> z_L_x_only = C_z_L->GetCompNonConst(0);
                C_z_L->Scal(1./rho);
                L1EPRestoNlp()->Px_L()->MultVector(-1., *tmpzL, 1., *tmp);

                SmartPtr<Vector> tmpzU = z_U->MakeNewCopy();
                auto C_z_U = dynamic_cast<CompoundVector*>(GetRawPtr(tmpzU));
                DBG_ASSERT(dynamic_cast<CompundVector*>(GetRawPtr(tmpzU)));
                SmartPtr<Vector> z_U_x_only = C_z_U->GetCompNonConst(0);
                C_z_U->Scal(1./rho);
                L1EPRestoNlp()->Px_U()->MultVector(1., *tmpzU, 1., *tmp);
            }
            else {
                L1EPRestoNlp()->Px_L()->MultVector(-1., *z_L, 1., *tmp);
                L1EPRestoNlp()->Px_U()->MultVector(1., *z_U, 1., *tmp);
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

    Number rho = L1EPRestoData().GetCurrentRho();
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

            if(L1EPRestoNlp()->l1_epr_inv_objective_type())
            {   // Take care of the (-z_L + z_U)/rho part.
                SmartPtr<Vector> tmpzL = z_L->MakeNewCopy();
                auto C_z_L = dynamic_cast<CompoundVector*>(GetRawPtr(tmpzL));
                DBG_ASSERT(dynamic_cast<CompundVector*>(GetRawPtr(tmpzL)));
                SmartPtr<Vector> z_L_x_only = C_z_L->GetCompNonConst(0);
                C_z_L->Scal(1./rho);
                L1EPRestoNlp()->Px_L()->MultVector(-1., *tmpzL, 1., *tmp);

                SmartPtr<Vector> tmpzU = z_U->MakeNewCopy();
                auto C_z_U = dynamic_cast<CompoundVector*>(GetRawPtr(tmpzU));
                DBG_ASSERT(dynamic_cast<CompundVector*>(GetRawPtr(tmpzU)));
                SmartPtr<Vector> z_U_x_only = C_z_U->GetCompNonConst(0);
                C_z_U->Scal(1./rho);
                L1EPRestoNlp()->Px_U()->MultVector(1., *tmpzU, 1., *tmp);
            }
            else {
                L1EPRestoNlp()->Px_L()->MultVector(-1., *z_L, 1., *tmp);
                L1EPRestoNlp()->Px_U()->MultVector(1., *z_U, 1., *tmp);
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

    Number rho = L1EPRestoData().GetCurrentRho();
    std::vector<Number> sdeps(1);
    sdeps[0] = rho;

    if( !curr_grad_lag_s_cache_l1_.GetCachedResult(result, deps, sdeps) )
    {
        if( !trial_grad_lag_s_cache_l1_.GetCachedResult(result, deps, sdeps) )
        {
            SmartPtr<Vector> tmp = y_d->MakeNew();
            L1EPRestoNlp()->Pd_U()->MultVector(1., *v_U, 0., *tmp);
            L1EPRestoNlp()->Pd_L()->MultVector(-1., *v_L, 1., *tmp);
            if(L1EPRestoNlp()->l1_epr_inv_objective_type())
            {
                tmp->Scal(1./rho);
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

    Number rho = L1EPRestoData().GetCurrentRho();
    std::vector<Number> sdeps(1);
    sdeps[0] = rho;

    if( !trial_grad_lag_s_cache_l1_.GetCachedResult(result, deps, sdeps) )
    {
        if( !curr_grad_lag_s_cache_l1_.GetCachedResult(result, deps, sdeps) )
        {
            SmartPtr<Vector> tmp = y_d->MakeNew();
            L1EPRestoNlp()->Pd_U()->MultVector(1., *v_U, 0., *tmp);
            L1EPRestoNlp()->Pd_L()->MultVector(-1., *v_L, 1., *tmp);
            if(L1EPRestoNlp()->l1_epr_inv_objective_type())
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
    Number rho = L1EPRestoData().GetCurrentRho();
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

        L1EPRestoNlp()->Px_L()->MultVector(kappa_d_l1_ * mu, *dampind_x_L, 1., *tmp);
        L1EPRestoNlp()->Px_U()->MultVector(-kappa_d_l1_ * mu, *dampind_x_U, 1., *tmp);

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
    Number rho = L1EPRestoData().GetCurrentRho();
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

        L1EPRestoNlp()->Pd_L()->MultVector(kappa_d_l1_ * mu, *dampind_s_L, 1., *tmp);
        L1EPRestoNlp()->Pd_U()->MultVector(-kappa_d_l1_ * mu, *dampind_s_U, 1., *tmp);

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
/*
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
*/
Vector& L1ExactPenaltyRestoCQ::Tmp_x_L_l1()
{
    if( !IsValid(tmp_x_L_l1_) )
    {
        tmp_x_L_l1_ = L1EPRestoNlp()->x_L()->MakeNew();
    }
    return *tmp_x_L_l1_;
}

Vector& L1ExactPenaltyRestoCQ::Tmp_x_U_l1()
{
    if( !IsValid(tmp_x_U_l1_) )
    {
        tmp_x_U_l1_ = L1EPRestoNlp()->x_U()->MakeNew();
    }
    return *tmp_x_U_l1_;
}

Vector& L1ExactPenaltyRestoCQ::Tmp_s_L_l1()
{
    if( !IsValid(tmp_s_L_l1_) )
    {
        tmp_s_L_l1_ = L1EPRestoNlp()->d_L()->MakeNew();
    }
    return *tmp_s_L_l1_;
}

Vector& L1ExactPenaltyRestoCQ::Tmp_s_U_l1()
{
    if( !IsValid(tmp_s_U_l1_) )
    {
        tmp_s_U_l1_ = L1EPRestoNlp()->d_U()->MakeNew();
    }
    return *tmp_s_U_l1_;
}

}