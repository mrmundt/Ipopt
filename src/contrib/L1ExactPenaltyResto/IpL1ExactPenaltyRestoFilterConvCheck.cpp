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
}


}