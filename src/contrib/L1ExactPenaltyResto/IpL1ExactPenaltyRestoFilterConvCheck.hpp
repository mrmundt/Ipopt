//
// Created by David on 1/3/2021.
//

#ifndef SRC_IPL1EXACTPENALTYRESTOFILTERCONVCHECK_HPP
#define SRC_IPL1EXACTPENALTYRESTOFILTERCONVCHECK_HPP

#include "IpOptErrorConvCheck.hpp"
#include "IpBacktrackingLSAcceptor.hpp"

namespace Ipopt
{

class L1ExactPenaltyRestoFilterConvCheck: public OptimalityErrorConvergenceCheck
{
public:
    L1ExactPenaltyRestoFilterConvCheck();

    ~L1ExactPenaltyRestoFilterConvCheck() override;

    bool InitializeImpl(
            const OptionsList& options,
            const std::string& prefix
            ) override;

    ConvergenceStatus CheckConvergence(
            bool call_intermediate_callback = true
            ) override;

    virtual void SetOirgLSAcceptor(
            const BacktrackingLSAcceptor& orig_ls_acceptor
            ) = 0;

    static void RegisterOptions(
            SmartPtr<RegisteredOptions> roptions
            );

private:

    L1ExactPenaltyRestoFilterConvCheck(
            const L1ExactPenaltyRestoFilterConvCheck&
            );

    void operator=(
            const L1ExactPenaltyRestoFilterConvCheck&
            );

    virtual ConvergenceStatus TestOrigProgress(
            Number orig_trial_barr,
            Number orit_trial_theta
            ) = 0;

    //Number kappa_resto_;

    Index maximum_iters_l1_;

    //Index maximum_resto_iters_;

    Number orig_constr_viol_tol_;

    bool first_resto_iter_;

    Index succesive_resto_iter_;
};

}
#endif //SRC_IPL1EXACTPENALTYRESTOFILTERCONVCHECK_HPP
