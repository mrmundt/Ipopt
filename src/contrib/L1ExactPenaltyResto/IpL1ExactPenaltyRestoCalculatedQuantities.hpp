//
// Created by David on 9/24/2020.
//

#ifndef SRC_IPL1EXACTPENALTYRESTOCALCULATEDQUANTITIES_HPP
#define SRC_IPL1EXACTPENALTYRESTOCALCULATEDQUANTITIES_HPP

#include "IpSmartPtr.hpp"
#include "IpCachedResults.hpp"

#include <string>

#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptNLP.hpp"
#include "IpL1ExactPenaltyRestoData.hpp"
#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"

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

    static void RegisterOptions(
            SmartPtr<RegisteredOptions> roptions
    );

    bool Initialize(
    const Journalist& jnlst,
    const OptionsList& options,
    const std::string& prefix
    );


    Number curr_f() override;
    Number trial_f() override;


    SmartPtr<const Vector> curr_grad_f() override;
    SmartPtr<const Vector> trial_grad_f() override;

    Number curr_barrier_obj() override;
    Number trial_barrier_obj() override;

    SmartPtr<const Vector> curr_grad_barrier_obj_x() override;
    SmartPtr<const Vector> curr_grad_barrier_obj_s() override;

    SmartPtr<const Vector> grad_kappa_times_damping_x() override;
    SmartPtr<const Vector> grad_kappa_times_damping_s() override;

    SmartPtr<const SymMatrix> curr_exact_hessian();

    SmartPtr<const Vector> curr_grad_lag_x() override;
    SmartPtr<const Vector> trial_grad_lag_x() override;

    SmartPtr<const Vector> curr_grad_lag_s() override;
    SmartPtr<const Vector> trial_grad_lag_s() override;

    SmartPtr<const Vector> curr_grad_lag_with_damping_x() override;
    SmartPtr<const Vector> curr_grad_lag_with_damping_s() override;


    /** Dual infeasibility in a given norm (at current iterate) */
    //virtual Number curr_dual_infeasibility(
    //        ENormType NormType
    //);
    /** Dual infeasibility in a given norm (at trial iterate) */
    //virtual Number trial_dual_infeasibility(
    //        ENormType NormType
    //);



    // curr_compl_x_U stays out


    SmartPtr<const Vector> curr_sigma_x() override;
    SmartPtr<const Vector> curr_sigma_s() override;


    Number ComputeRhoTrial();
private:
    /**@name Default Compiler Generated Methods
     *
     */
    L1ExactPenaltyRestoCQ();

    L1ExactPenaltyRestoCQ(const L1ExactPenaltyRestoCQ&);

    void operator=(const L1ExactPenaltyRestoCQ&);

    enum RhoUpdateKind{
        QUAD=0,
        QUADNOSIGNA,
        LINEAR,
        CONST
    };

    RhoUpdateKind l1exactpenalty_rho_type_;

    SmartPtr<IpoptData> ip_data_l1_;
    SmartPtr<IpoptNLP> ip_nlp_l1_;
    SmartPtr<L1ExactPenaltyRestoCQ> l1epr_cq_;

    L1ExactPenaltyRestoData& L1EPRestoData()
    {
        L1ExactPenaltyRestoData* l1epresto_data = static_cast<L1ExactPenaltyRestoData*>(&(ip_data_l1_->AdditionalData()));
        DBG_ASSERT(dynamic_cast<L1ExactPenaltyRestoData*>(&(ip_data_l1->AdditionalData())));
        return l1epresto_data;
    }

    SmartPtr<IpL1ExactPenaltyRestoIpoptNLP> L1EPRestoNlp()
    {
        SmartPtr<const IpL1ExactPenaltyRestoIpoptNLP> l1epr_nlp = static_cast<const IpL1ExactPenaltyRestoIpoptNLP*>(GetRawPtr(ip_nlp_l1_));
        DBG_ASSERT(dynamic_cast<const IpL1ExactPenaltyRestoIpoptNLP*>(GetRawPtr(ip_nlp_l1_)));
        return l1epr_nlp;
    }
    SmartPtr<Vector> tmp_x_l1_;
    SmartPtr<Vector> tmp_s_l1_;
    SmartPtr<Vector> tmp_c_l1_;
    SmartPtr<Vector> tmp_d_l1_;
    SmartPtr<Vector> tmp_x_L_l1_;
    SmartPtr<Vector> tmp_x_U_l1_;
    SmartPtr<Vector> tmp_s_L_l1_;
    SmartPtr<Vector> tmp_s_U_l1_;

    /** Accessor methods for the temporary vectors */
    Vector& Tmp_x_l1();
    Vector& Tmp_s_l1();
    Vector& Tmp_c_l1();
    Vector& Tmp_d_l1();
    Vector& Tmp_x_L_l1();
    Vector& Tmp_x_U_l1();
    Vector& Tmp_s_L_l1();
    Vector& Tmp_s_U_l1();


    Number CalcBarrierTermL1(
            Number mu,
            const Vector& slack_x_L,
            const Vector& slack_x_U,
            const Vector& slack_s_L,
            const Vector& slack_s_U
    );

    void ComputeDampingIndicatorsL1(
            SmartPtr<const Vector>& dampind_x_L,
            SmartPtr<const Vector>& dampind_x_U,
            SmartPtr<const Vector>& dampind_s_L,
            SmartPtr<const Vector>& dampind_s_U
    );

    SmartPtr<Vector> CalcComplL1(
            const Vector& slack,
            const Vector& mult
    );

    CachedResults<Number> curr_f_cache_l1_;
    CachedResults<Number> trial_f_cache_l1_;
    CachedResults<const Vector> curr_grad_f_cache_l1_;
    CachedResults<const Vector> trial_grad_f_cache_l1_;
    CachedResults<Number> curr_barrier_obj_cache_l1_;
    CachedResults<Number> trial_barrier_obj_cache_l1_;

    CachedResults<const Vector> curr_grad_barrier_obj_x_cache_l1_;
    CachedResults<const Vector> curr_grad_barrier_obj_s_cache_l1_;

    CachedResults<const Vector> grad_kappa_times_damping_x_cache_l1_;
    CachedResults<const Vector> grad_kappa_times_damping_s_cache_l1_;


    CachedResults<const Vector> dampind_x_L_cache_l1_;
    CachedResults<const Vector> dampind_x_U_cache_l1_;
    CachedResults<const Vector> dampind_s_L_cache_l1_;
    CachedResults<const Vector> dampind_s_U_cache_l1_;


    SmartPtr<Vector> dampind_x_L_l1_;
    SmartPtr<Vector> dampind_x_U_l1_;
    SmartPtr<Vector> dampind_s_L_l1_;
    SmartPtr<Vector> dampind_s_U_l1_;

    Number kappa_d_l1_;

    CachedResults<const SymMatrix> curr_exact_hessian_cache_l1_;
    CachedResults<const Vector> curr_grad_lag_x_cache_l1_;
    CachedResults<const Vector> trial_grad_lag_x_cache_l1_;
    CachedResults<const Vector> curr_grad_lag_s_cache_l1_;
    CachedResults<const Vector> trial_grad_lag_s_cache_l1_;

    CachedResults<const Vector> curr_grad_lag_with_damping_x_cache_l1_;
    CachedResults<const Vector> curr_grad_lag_with_damping_s_cache_l1_;

    CachedResults<const Vector> curr_compl_x_L_cache_l1_;
    CachedResults<const Vector> trial_compl_x_L_cache_l1_;

    CachedResults<const Vector> curr_compl_x_U_cache_l1_;
    CachedResults<const Vector> trial_compl_x_U_cache_l1_;

    CachedResults<const Vector> curr_compl_s_L_cache_l1_;
    CachedResults<const Vector> curr_compl_s_U_cache_l1_;

    CachedResults<const Vector> trial_compl_s_L_cache_l1_;
    CachedResults<const Vector> trial_compl_s_U_cache_l1_;

    CachedResults<const Vector> curr_relaxed_compl_x_L_cached_l1_;
    CachedResults<const Vector> curr_relaxed_compl_x_U_cached_l1_;

    CachedResults<const Vector> curr_relaxed_compl_s_L_cached_l1_;
    CachedResults<const Vector> curr_relaxed_compl_s_U_cached_l1_;

    CachedResults<const Vector> curr_complementarity_cache_l1_;
    CachedResults<const Vector> trial_complementarity_cache_l1_;

    CachedResults<const Vector> curr_sigma_x_cache_l1_;
    CachedResults<const Vector> curr_sigma_s_cache_l1_;


};
}

#endif //SRC_IPL1EXACTPENALTYRESTOCALCULATEDQUANTITIES_HPP
