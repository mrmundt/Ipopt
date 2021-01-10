//
// Created by David on 1/9/2021.
//

#ifndef SRC_IPL1EXACTPENALTYRHO_HPP
#define SRC_IPL1EXACTPENALTYRHO_HPP

#include "IpAlgStrategy.hpp"

namespace Ipopt
{
class L1ExactPenaltyRho : public AlgorithmStrategyObject
{
public:
    L1ExactPenaltyRho();

    ~L1ExactPenaltyRho();

    bool InitializeImpl(
            const OptionsList& options,
            const std::string& prefix
            ) override;

    Number ComputeRhoTrial();

private:
    L1ExactPenaltyRho(
            const L1ExactPenaltyRho&
            );

    void operator=(
            const L1ExactPenaltyRho&
            );

    CachedResults<Number> trial_rho_cache_;


};
}



#endif //SRC_IPL1EXACTPENALTYRHO_HPP
