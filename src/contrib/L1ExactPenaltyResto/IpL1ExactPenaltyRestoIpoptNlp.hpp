//
// Created by David on 10/1/2020.
//

#ifndef SRC_IPL1EXACTPENALTYRESTOIPOPTNLP_HPP
#define SRC_IPL1EXACTPENALTYRESTOIPOPTNLP_HPP

#include "IpRestoIpoptNLP.hpp"

namespace Ipopt
{


class L1ExactPenaltyRestoIpoptNLP : public RestoIpoptNLP
{
public:
    L1ExactPenaltyRestoIpoptNLP(
            IpoptNLP&   orig_ip_nlp,
            IpoptData&  orig_ip_data,
            IpoptCalculatedQuantities&  orig_ip_cq,
            const SmartPtr<IpoptData>& l1_ip_data
    );

    ~L1ExactPenaltyRestoIpoptNLP();

    bool Initialize(
        const Journalist& jnlst,
        const OptionsList& options,
        const std::string& prefix
    ) override;

    static void RegisterOptions(
        SmartPtr<RegisteredOptions> roptions
    );

    enum IpL1ExactPenaltyObjectiveType{
        CONSTRAINT=0,
        OBJECTIVE_INV
    };

    bool l1_epr_inv_objective_type() const;
    //l1_epr_suff_feasib_update_

    Number f(const Vector &x) override;
    Number f(const Vector &x, Number rho) override;

    SmartPtr<const Vector> grad_f(const Vector &x) override;
    SmartPtr<const Vector> grad_f(const Vector &x, Number rho) override;
    // so I figured these functions either depend on mu or rho but not both.
    SmartPtr<const Vector> c(const Vector &x) override;

    SmartPtr<const SymMatrix> h(
            const Vector &x,
            Number obj_factor,
            const Vector &yc,
            const Vector &yd) override;

    SmartPtr<const SymMatrix> h(
            const Vector &x,
            Number obj_factor,
            const Vector &yc,
            const Vector &yd,
            Number rho) override;


    SmartPtr<const SymMatrix> uninitialized_h() override;
    bool objective_depends_on_mu_rho() const
    {
        return true;
    }

    Number Rho() const override;


private:
    IpL1ExactPenaltyObjectiveType l1_epr_objective_type_{CONSTRAINT};
    SmartPtr<IpoptData> l1_ip_data_;
    SmartPtr<Vector> l1_diag_vec_dummy_;
    SmartPtr<DiagMatrix> l1_diag_mat_dum_;
    HessianApproximationType hessian_approximation_l1_{EXACT};
public:
    SmartPtr<const DiagMatrix> getL1DiagMatDummy();

public:
    IpL1ExactPenaltyObjectiveType getL1EprObjectiveType() const {
        return l1_epr_objective_type_;
    }
};
}



#endif //SRC_IPL1EXACTPENALTYRESTOIPOPTNLP_HPP
