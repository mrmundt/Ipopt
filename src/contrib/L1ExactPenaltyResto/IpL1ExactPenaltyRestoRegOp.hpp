//
// Created by David on 12/15/2020.
//

#ifndef SRC_IPL1EXACTPENALTYRESTOREGOP_HPP
#define SRC_IPL1EXACTPENALTYRESTOREGOP_HPP

#include "IpSmartPtr.hpp"

namespace Ipopt
{
    class RegisteredOptions;

    void RegisterOptions_L1ExactPenaltyResto(
        const SmartPtr<RegisteredOptions>& roptions
    );
}

#endif //SRC_IPL1EXACTPENALTYRESTOREGOP_HPP
