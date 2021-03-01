//
// Created by dav0 on 2/20/21.
//

#ifndef IPOPT_L1EPR_IPL1EXACTPENALTYLEASTSQUAREMULTS_HPP
#define IPOPT_L1EPR_IPL1EXACTPENALTYLEASTSQUAREMULTS_HPP

#include "IpEqMultCalculator.hpp"
#include "IpAugRestoSystemSolver.hpp"
#include "IpL1ExactPenaltyRestoData.hpp"

namespace Ipopt
{
class L1ExactPenaltyLeastSquareMults : public EqMultiplierCalculator {
public:
    explicit L1ExactPenaltyLeastSquareMults(AugSystemSolver& augSystemSolver);
    ~L1ExactPenaltyLeastSquareMults() override = default;
    bool InitializeImpl(
            const OptionsList& options,
            const std::string& prefix
            ) override;
    bool CalculateMultipliers(
            Vector& y_c,
            Vector& y_d
            ) override;
    L1ExactPenaltyLeastSquareMults() = delete;
    L1ExactPenaltyLeastSquareMults(const L1ExactPenaltyLeastSquareMults&) = delete;
    L1ExactPenaltyLeastSquareMults& operator=(const L1ExactPenaltyLeastSquareMults&) = delete;

private:
    SmartPtr<AugSystemSolver> augsyssolver_;
    bool is_rho_inv_obj_{false};
    L1ExactPenaltyRestoData* l1EprAddData_{nullptr};
};

}


#endif //IPOPT_L1EPR_IPL1EXACTPENALTYLEASTSQUAREMULTS_HPP
