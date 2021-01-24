//
// Created by David on 12/15/2020.
//
#include "IpAlgorithmRegOp.hpp"
#include "IpRegOptions.hpp"

#include "IpL1ExactPenaltyRestoRegOp.hpp"
#include "IpL1ExactPenaltyRestoData.hpp"
#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"
#include "IpL1ExactPenaltyRhoUpdater.hpp"
#include "IpL1ExactPenaltyRestoCalculatedQuantities.hpp"

namespace Ipopt
{
    void RegisterOptions_L1ExactPenaltyResto(
        const SmartPtr<RegisteredOptions>& roptions
    )
{
    roptions->SetRegisteringCategory("l1 Exact Penalty");
    L1ExactPenaltyRestoData::RegisterOptions(roptions);
    L1ExactPenaltyRestoCQ::RegisterOptions(roptions);
    L1ExactPenaltyRestoIpoptNLP::RegisterOptions(roptions);
    L1ExactPenaltyRhoUpdater::RegisterOptions(roptions);

}
}
