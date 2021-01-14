//
// Created by David on 11/23/2020.
//

#include "IpL1ExactPenaltyRestoData.hpp"

namespace Ipopt
{
L1ExactPenaltyRestoData::L1ExactPenaltyRestoData()
{ }

L1ExactPenaltyRestoData::~L1ExactPenaltyRestoData() noexcept {}


bool L1ExactPenaltyRestoData::Initialize(const Journalist &jnlst,
                                         const OptionsList &options,
                                         const std::string &prefix) {
    return true;
}

bool L1ExactPenaltyRestoData::InitializeDataStructures() {
    rho_initialized_ = true;
    rho_ = 1.00;
    rho_trial_ = 1.00;
    rho_has_changed_ = false;
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