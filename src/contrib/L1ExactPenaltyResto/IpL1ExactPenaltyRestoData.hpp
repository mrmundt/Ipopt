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
    L1ExactPenaltyRestoData();

    virtual ~L1ExactPenaltyRestoData() override;

    bool Initialize(
            const Journalist& jnlst,
            const OptionsList& options,
            const std::string& prefix
            ) override;

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
    Number rho_;
    Number rho_trial_;
    bool rho_initialized_ = true;
    bool rho_has_changed_ = true;

    Number SetRho(Number rho) {
        rho_ = rho;
        rho_initialized_ = true;
        rho_has_changed_ = true;
    }

};
}


#endif //SRC_IPL1EXACTPENALTYRESTODATA_HPP
