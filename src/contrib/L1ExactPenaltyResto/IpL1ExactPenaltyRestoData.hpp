//
// Created by David on 11/23/2020.
//

#ifndef SRC_IPL1EXACTPENALTYRESTODATA_HPP
#define SRC_IPL1EXACTPENALTYRESTODATA_HPP

#include "IpIpoptData.hpp"
#include "IpOptionsList.hpp"

namespace Ipopt
{
class L1ExactPenaltyRestoData: public IpoptAdditionalData
{
public:
    L1ExactPenaltyRestoData() = default;

    ~L1ExactPenaltyRestoData() override = default;
    L1ExactPenaltyRestoData(const L1ExactPenaltyRestoData&) = delete;
    L1ExactPenaltyRestoData& operator=(const L1ExactPenaltyRestoData&) = delete;

    bool Initialize(
            const Journalist& jnlst,
            const OptionsList& options,
            const std::string& prefix
            ) override;

    static void RegisterOptions(
            SmartPtr<RegisteredOptions> roptions
            );

    bool InitializeDataStructures() override;
    void AcceptTrialPoint() override;

    Number GetCurrentRho() const {
        return rho_;
    }

    void SetRhoTrial(Number rho){
        rho_trial_ = rho;
    }

    void SetRhoStatus(bool stat) {
        rho_has_changed_ = stat;
    }

    bool GetRhoStatus() const
    {
        return rho_has_changed_;
    }

    void AcceptRhoTrial();

private:
    Number rho_init_{1.};
    Number rho_{1.};
    Number rho_trial_{1.};
    bool rho_initialized_ = true;
    bool rho_has_changed_ = true;

    void SetRho(Number rho) {
        rho_ = rho;
        rho_initialized_ = true;
        rho_has_changed_ = true;
    }

};
}


#endif //SRC_IPL1EXACTPENALTYRESTODATA_HPP
