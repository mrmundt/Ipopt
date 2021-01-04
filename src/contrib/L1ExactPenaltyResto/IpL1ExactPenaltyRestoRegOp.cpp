//
// Created by David on 12/15/2020.
//
#include "IpAlgorithmRegOp.hpp"
#include "IpRegOptions.hpp"

#include "IpL1ExactPenaltyRestoRegOp.hpp"
#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"
#include "IpL1ExactPenaltyRestoCalculatedQuantities.hpp"

namespace Ipopt
{
    void RegisterOptions_L1ExactPenaltyResto(
        const SmartPtr<RegisteredOptions>& roptions;
    )
{
    roptions->SetRegisteringCategory("l1 Exact Penalty");
    L1ExactPenaltyRestoCQ::RegisterOptions(roptions);
    IpL1ExactPenaltyRestoIpoptNLP::RegisterOptions(roptions);

}
}
