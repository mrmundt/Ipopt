//
// Created by David on 9/24/2020.
//

#ifndef SRC_IPL1EXACTPENALTYRESTOCALCULATEDQUANTITIES_HPP
#define SRC_IPL1EXACTPENALTYRESTOCALCULATEDQUANTITIES_HPP

#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptNLP.hpp"
#include "IpL1ExactPenaltyRestoData.hpp"
#include "IpCachedResults.hpp"

namespace Ipopt
{
/** Class for the calculated quantities for the l1-exact penalty formulation
 *
 */
class L1ExactPenaltyRestoCQ: public IpoptCalculatedQuantities
{
public:
    /**@name Constuctors/Destructors */
    //@{
    /** Default Constructor */
    L1ExactPenaltyRestoCQ(
            const SmartPtr<IpoptNLP>& ip_nlp,
            const SmartPtr<IpoptData>& ip_data
            );

    /** Destructor */
    virtual ~L1ExactPenaltyRestoCQ();
    //@}

    Number curr_f() override;

    Number trial_f() override;


    SmartPtr<const Vector> curr_grad_f() override;

    SmartPtr<const Vector> trial_grad_f() override;

    Number curr_barrier_obj() override;
    Number trial_barrier_obj() override;

    SmartPtr<const Vector> curr_grad_barrier_obj_x() override;
    SmartPtr<const Vector> curr_grad_barrier_obj_s() override;

    SmartPtr<const SymMatrix> curr_exact_hessian();

    SmartPtr<const Vector> curr_grad_lag_x();
    SmartPtr<const Vector> trial_grad_lag_x();

    SmartPtr<const Vector> curr_grad_lag_with_damping_x();
    SmartPtr<const Vector> curr_grad_lag_with_damping_s();
    // curr_compl_x_U stays out
    Number curr_barrier_error();

    SmartPtr<const Vector> curr_sigma_x();

    Number curr_gradBarrTDelta();

    Number ComputeRhoTrial();
private:
    /**@name Default Compiler Generated Methods
     *
     */
    L1ExactPenaltyRestoCQ();

    L1ExactPenaltyRestoCQ(const L1ExactPenaltyRestoCQ&);

    void operator=(const L1ExactPenaltyRestoCQ&);

    SmartPtr<IpoptData> ip_data_l1_;
    SmartPtr<IpoptNLP> ip_nlp_l1_;

    L1ExactPenaltyRestoData& L1EPRestoData()
    {
        L1ExactPenaltyRestoData& l1epresto_data = static_cast<L1ExactPenaltyRestoData&>(ip_data_l1_->AdditionalData());
        DBG_ASSERT(dynamic_cast<L1ExactPenaltyRestoData*>(&ip_data_l1->AdditionalData()));
        return l1epresto_data;
    }


    Number CalcBarrierTerm(
            Number mu,
            const Vector& slack_x_L,
            const Vector& slack_x_U,
            const Vector& slack_s_L,
            const Vector& slack_s_U
    );

    void ComputeDampingIndicators(
            SmartPtr<const Vector>& dampind_x_L,
            SmartPtr<const Vector>& dampind_x_U,
            SmartPtr<const Vector>& dampind_s_L,
            SmartPtr<const Vector>& dampind_s_U
    );

    void ScaledDamping_x(SmartPtr<Vector> &d1,
                         SmartPtr<Vector> &d2,
                         Number rho,
                         bool scale_rho);

    CachedResults<Number> curr_f_cache_l1_;
    CachedResults<Number> trial_f_cache_l1_;


};
}

#endif //SRC_IPL1EXACTPENALTYRESTOCALCULATEDQUANTITIES_HPP
