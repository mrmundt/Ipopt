//
// Created by David on 10/12/2020.
//

#include "IpL1ExactPenaltyRestoCalculatedQuantities.hpp"

namespace Ipopt
{
//#if IPOPT_VERBOSITY > 0
    static const Index dbg_verbosity = 0;
//#endif
L1ExactPenaltyRestoCQ::L1ExactPenaltyRestoCQ(
        const SmartPtr<IpoptNLP>& ip_nlp,
        const SmartPtr<IpoptData>& ip_data
        )
        : IpoptCalculatedQuantities(ip_nlp, ip_data),
        ip_nlp_l1_(ip_nlp),
        ip_data_l1_(ip_data),
        curr_f_cache_l1_(1)
{
    DBG_ASSERT(IsValid(ip_data_));
}
L1ExactPenaltyRestoCQ::~L1ExactPenaltyRestoCQ() noexcept
{ }

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
    Number rho = l1data.CurrentRho();
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

}

}