//
// Created by David on 1/9/2021.
//

#ifndef SRC_IPL1EXACTPENALTYRHOUPDATER_HPP
#define SRC_IPL1EXACTPENALTYRHOUPDATER_HPP



#include "IpAlgStrategy.hpp"
#include "IpSumSymMatrix.hpp"
#include "IpDiagMatrix.hpp"
#include "IpL1ExactPenaltyRestoData.hpp"


namespace Ipopt
{
class L1ExactPenaltyRhoUpdater : public AlgorithmStrategyObject
{
public:
    L1ExactPenaltyRhoUpdater();

    ~L1ExactPenaltyRhoUpdater();

    bool InitializeImpl(
            const OptionsList& options,
            const std::string& prefix
            ) override;



    void UpdateRhoTrial();
    void UpdateRhoAction();
    enum RhoUpdateKind{
        QUAD=0,
        QUADNOSIGNA,
        LINEAR,
        CONST
    };

private:
    L1ExactPenaltyRhoUpdater(
            const L1ExactPenaltyRhoUpdater&
            );

    void operator=(
            const L1ExactPenaltyRhoUpdater&
            );

    CachedResults<Number> trial_rho_cache_;

    RhoUpdateKind l1_epr_update_kind_;
    Number l1_epr_rho0_;
    Number l1_epr_epsi_;
    Number l1_epr_max_rho;

    bool l1_epr_suff_feasib_update_;
    bool l1_epr_has_changed_;
    SmartPtr<DiagMatrixSpace> Sigma_x_space_;
    SmartPtr<SumSymMatrixSpace> Hx_p_Sigma_x_space_;

    SmartPtr<const SymMatrixSpace> original_h_space_;

    SmartPtr<SumSymMatrix> Hx_Sigma_x_;
    SmartPtr<DiagMatrix> Sigma_x_;

    void SetMatrixSpaces();

    SumSymMatrix& Tmp_Wx_Sigma_x();
    DiagMatrix& Tmp_Sigma_x();
    Number ComputeRhoTrial();
    L1ExactPenaltyRestoData& L1EPRAddData();


};
}



#endif //SRC_IPL1EXACTPENALTYRHOUPDATER_HPP
