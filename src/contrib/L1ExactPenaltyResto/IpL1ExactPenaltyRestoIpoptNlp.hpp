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
            IpoptCalculatedQuantities&  orig_ip_cq
    );

    ~L1ExactPenaltyRestoIpoptNLP() noexcept;

    virtual bool Initialize(
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

    bool l1exactpenalty_inv_objective_type() const;

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




private:
    IpL1ExactPenaltyObjectiveType l1_exact_penalty_objective_type_;
};
}



#endif //SRC_IPL1EXACTPENALTYRESTOIPOPTNLP_HPP
