//
// Created by David on 1/9/2021.
//

#ifndef SRC_IPL1EXACTPENALTYRHOUPDATER_HPP
#define SRC_IPL1EXACTPENALTYRHOUPDATER_HPP



#include "IpAlgStrategy.hpp"
#include "IpSumSymMatrix.hpp"
#include "IpDiagMatrix.hpp"
#include "IpL1ExactPenaltyRestoData.hpp"
#include "IpBacktrackingLSAcceptor.hpp"

namespace Ipopt
{
class L1ExactPenaltyRhoUpdater : public AlgorithmStrategyObject
{
public:
    explicit L1ExactPenaltyRhoUpdater(
            const SmartPtr<BacktrackingLSAcceptor>& resto_bls_acceptor
            );

    ~L1ExactPenaltyRhoUpdater() override = default;

    bool InitializeImpl(
            const OptionsList& options,
            const std::string& prefix
            ) override;

    static void RegisterOptions(
            SmartPtr<RegisteredOptions> roptions
            );

    bool UpdateRhoTrial();
    void UpdateRhoAction();
    enum RhoUpdateKind{
        QUAD=0,
        QUADNOSIGNA,
        LINEAR,
        CONST
    };

    L1ExactPenaltyRhoUpdater() = delete;
    L1ExactPenaltyRhoUpdater(
            const L1ExactPenaltyRhoUpdater&
    ) = delete;

    void operator=(
            const L1ExactPenaltyRhoUpdater&
    ) = delete;

private:

    SmartPtr<BacktrackingLSAcceptor> ip_bls_acceptor_;

    CachedResults<Number> trial_rho_cache_;

    RhoUpdateKind l1_epr_update_kind_;
    Number l1_epr_epsi_{1.};
    Number l1_epr_max_rho{1e+07};

    bool l1_epr_suff_feasib_update_{false};
    bool l1_epr_has_changed_{false};
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
    std::string resto_lsacceptor_option_;


};
}



#endif //SRC_IPL1EXACTPENALTYRHOUPDATER_HPP
