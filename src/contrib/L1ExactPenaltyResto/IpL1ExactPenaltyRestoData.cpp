//
// Created by David on 11/23/2020.
//

#include "IpL1ExactPenaltyRestoData.hpp"

namespace Ipopt
{
L1ExactPenaltyRestoData::L1ExactPenaltyRestoData()
= default;

L1ExactPenaltyRestoData::~L1ExactPenaltyRestoData() noexcept = default;

void L1ExactPenaltyRestoData::RegisterOptions(
        SmartPtr<RegisteredOptions> roptions)
{
    roptions->AddLowerBoundedNumberOption(
            "l1_init_penalty",
            "Initial value of rho.",
            0.0, true, 1e+03,
            "Initial value of rho for the l1-exact penalty formulation");

}


bool L1ExactPenaltyRestoData::Initialize(const Journalist &jnlst,
                                         const OptionsList &options,
                                         const std::string &prefix) {

    rho_init_ = 1.;
    options.GetNumericValue("l1_init_penalty", rho_init_, prefix);
    return true;
}

bool L1ExactPenaltyRestoData::InitializeDataStructures() {
    rho_initialized_ = true;
    rho_ = rho_init_;
    rho_trial_ = rho_init_;
    rho_has_changed_ = false;
    return true;
}

void L1ExactPenaltyRestoData::AcceptTrialPoint()
{

}

void L1ExactPenaltyRestoData::AcceptRhoTrial()
{
    rho_ = rho_trial_;
    rho_trial_ = 0.;
}


}