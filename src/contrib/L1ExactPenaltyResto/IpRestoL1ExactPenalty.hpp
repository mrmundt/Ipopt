//
// Created by David on 1/14/2021.
//

#ifndef SRC_IPRESTOL1EXACTPENALTY_HPP
#define SRC_IPRESTOL1EXACTPENALTY_HPP

#include "IpRestoPhase.hpp"
#include "IpL1IpoptAlg.hpp"
#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"
#include "IpEqMultCalculator.hpp"

namespace Ipopt {
class L1ExactPenaltyRestorationPhase: public RestorationPhase
{
public:
    L1ExactPenaltyRestorationPhase(
            L1IpoptAlg&   resto_alg,
            const SmartPtr<EqMultiplierCalculator>& eq_mult_calculator
            );

    ~L1ExactPenaltyRestorationPhase() override;

    bool InitializeImpl(const OptionsList &options, const std::string &prefix) override;

    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);

protected:
    bool PerformRestoration() override;

private:
    L1ExactPenaltyRestorationPhase();

    L1ExactPenaltyRestorationPhase(const L1ExactPenaltyRestorationPhase&);

    void operator=(const L1ExactPenaltyRestorationPhase&);

    SmartPtr<L1IpoptAlg> resto_alg_;
    SmartPtr<EqMultiplierCalculator> eq_mult_calculator_;
    SmartPtr<OptionsList> resto_options_;
    Number constr_mult_reset_threshold_;
    Number bound_mult_reset_threshold_;


    Number constr_viol_tol_;
    Number resto_failure_feasibility_threshold_;
    Index count_restorations_;

    void ComputeBoundMultiplierStep(
            Vector&       delta_z,
            const Vector& curr_z,
            const Vector& curr_slack,
            const Vector& trial_slack
            );

};
}

#endif //SRC_IPRESTOL1EXACTPENALTY_HPP
