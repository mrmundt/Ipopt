//
// Created by David on 1/9/2021.
//

#ifndef SRC_IPL1EXACTPENALTYRHO_HPP
#define SRC_IPL1EXACTPENALTYRHO_HPP

#include "IpAlgStrategy.hpp"
#include "IpSumSymMatrix.hpp"
#include "IpDiagMatrix.hpp"

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

    enum RhoUpdateKind{
        QUAD=0,
        QUADNOSIGNA,
        LINEAR,
        CONST
    };

    RhoUpdateKind l1_epr_update_kind_;
    SmartPtr<DiagMatrixSpace> Sigma_x_space_;
    SmartPtr<SumSymMatrixSpace> Hx_p_Sigma_x_space_;

    SmartPtr<const SymMatrixSpace> original_h_space_;

    SmartPtr<SumSymMatrix> Hx_Sigma_x_;
    SmartPtr<DiagMatrix> Sigma_x_;

    void SetMatrixSpaces();

    SumSymMatrix& Tmp_Wx_Sigma_x();
    DiagMatrix& Tmp_Sigma_x();

};
}



#endif //SRC_IPL1EXACTPENALTYRHO_HPP
