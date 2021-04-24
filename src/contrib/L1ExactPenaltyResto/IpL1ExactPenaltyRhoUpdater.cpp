//
// Created by David on 1/9/2021.
//

#include "IpL1ExactPenaltyRhoUpdater.hpp"
#include "IpCompoundVector.hpp"
#include "IpCompoundSymMatrix.hpp"
#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"
#include "IpFilterLSAcceptor.hpp"


namespace Ipopt
{

L1ExactPenaltyRhoUpdater::L1ExactPenaltyRhoUpdater(
        const SmartPtr<BacktrackingLSAcceptor> &resto_bls_acceptor) :
        ip_bls_acceptor_(resto_bls_acceptor),
        trial_rho_cache_(2),
        l1_epr_update_kind_(LINEAR)
{
    DBG_START_METH("L1ExactPenaltyRhoUpdater::L1ExactPenaltyRhoUpdater()", dbg_verbosity);
}

void L1ExactPenaltyRhoUpdater::RegisterOptions(
        SmartPtr<RegisteredOptions> roptions) {
    roptions->SetRegisteringCategory("l1 Exact Penalty");
    roptions->AddStringOption4(
            "l1_penalty_type",
            "Type of update for the penalty parameter",
            "linear_model",
            "quadratic_model", "check the quadratic model",
            "quadratic_model_no_sigma", "quadratic model without the barrier",
            "linear_model", "use a linear model for predicted reduction",
            "fixed", "fixed rho",
            "Type of update for the penalty."
    );
    roptions->AddBoundedNumberOption(
            "l1_epsilon",
            "eps of (1 - eps)||c|| - (p+ + n+)",
            0.0, true, 1.0, true, 0.1,
            "Determines the aggressiveness of the denominator for rho");
    roptions->AddLowerBoundedNumberOption(
            "l1_penalty_feas_max",
            "Max value of the rho update by factor.",
            0.0, true, 1e+07,
            "Set this value to the maximum allowed rho if the Factor"
            "update strategy is used");
}


bool L1ExactPenaltyRhoUpdater::InitializeImpl(const OptionsList &options,
                                              const std::string &prefix) {
    SetMatrixSpaces();
    options.GetStringValue("line_search_method", resto_lsacceptor_option_, "resto." + prefix);
    Index l1rhotype;
    options.GetEnumValue("l1_penalty_type", l1rhotype, prefix);
    l1_epr_update_kind_ = RhoUpdateKind(l1rhotype);
    options.GetNumericValue("l1_epsilon", l1_epr_epsi_, prefix);
    options.GetNumericValue("l1_penalty_feas_max", l1_epr_max_rho, prefix);
    return true;
}


Number L1ExactPenaltyRhoUpdater::ComputeRhoTrial()
{

    if (l1_epr_update_kind_ == CONST)
    {
        return L1EPRAddData().GetCurrentRho();
    } // if fixed

    // rho
    Number rho_trial;
    SmartPtr<const Vector> x_resto = IpData().curr()->x();
    SmartPtr<const Vector> s_resto = IpData().curr()->s();
    SmartPtr<const Vector> dx_R = IpData().delta()->x();
    SmartPtr<const Vector> ds = IpData().delta()->s();

    std::vector<const TaggedObject*> tdeps(4);
    tdeps[0] = GetRawPtr(x_resto);
    tdeps[1] = GetRawPtr(s_resto);
    tdeps[2] = GetRawPtr(dx_R);
    tdeps[3] = GetRawPtr(ds);
    Number mu = IpData().curr_mu();
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;
    trial_rho_cache_.GetCachedResult(rho_trial, tdeps, sdeps);

    SmartPtr<const CompoundVector> xC = dynamic_cast<const CompoundVector*>(GetRawPtr(x_resto));

    // fetch grad_f_barrier
    SmartPtr<const Vector> grad_barrier_xR = IpCq().curr_grad_barrier_obj_x();
    SmartPtr<const CompoundVector> grad_barrier_xRC = dynamic_cast<const CompoundVector*>(GetRawPtr(grad_barrier_xR));
    SmartPtr<const Vector> grad_phi_x = grad_barrier_xRC->GetComp(0);
    SmartPtr<const Vector> grad_phi_s = IpCq().curr_grad_barrier_obj_s();
    // fetch dx
    SmartPtr<const CompoundVector> dC = dynamic_cast<const CompoundVector*>(GetRawPtr(dx_R));
    SmartPtr<const Vector> dx = dC->GetComp(0);
    ////
    // fetch p and n
    SmartPtr<Vector> nplusc = xC->GetComp(1)->MakeNewCopy();
    SmartPtr<Vector> pplusc = xC->GetComp(2)->MakeNewCopy();
    SmartPtr<Vector> nplusd = xC->GetComp(3)->MakeNewCopy();
    SmartPtr<Vector> pplusd = xC->GetComp(4)->MakeNewCopy();
    // debug
    //SmartPtr<Vector> dummy = xC->GetComp(1)->MakeNew();
    //dummy->Copy(*nplusc);
    //dummy->ElementWiseMultiply(*pplusc);
    // curr c and d
    SmartPtr<const Vector> c_R = IpCq().curr_c();
    SmartPtr<const Vector> d_minus_s_R = IpCq().curr_d_minus_s();
    SmartPtr<const CompoundVector> cC = dynamic_cast<const CompoundVector*>(GetRawPtr(c_R));
    SmartPtr<const CompoundVector> dmsC = dynamic_cast<const CompoundVector*>(GetRawPtr(d_minus_s_R));

    SmartPtr<Vector> c = cC->GetComp(0)->MakeNewCopy();
    SmartPtr<Vector> dms = dmsC->GetComp(0)->MakeNewCopy();

    c->AddTwoVectors(-1., *nplusc, 1.,  *pplusc, 1.);
    dms->AddTwoVectors(-1., *nplusd, 1., *pplusd, 1.);


    // ck
    Number ck = IpCq().CalcNormOfType(NORM_1, *c, *dms);

    // calculate p+ and n+
    nplusc->AddOneVector(1., *(dC->GetComp(1)), 1.);
    pplusc->AddOneVector(1., *(dC->GetComp(2)), 1.);
    nplusd->AddOneVector(1., *(dC->GetComp(3)), 1.);
    pplusd->AddOneVector(1., *(dC->GetComp(4)), 1.);
    //DBG_PRINT_VECTOR(0, "dummy", *dummy);
    pplusc->AddOneVector(+1., *nplusc, 1.0);
    pplusd->AddOneVector(+1., *nplusd, 1.0);
    Number ppnp = pplusc->Sum() + pplusd->Sum();
    ////
    Number dxWxdx = 0.;
    Number dsWsdx = 0.;
    Number sigW = 0.;
    //THROW_EXCEPTION(DAVS, "fail\n");
    //ASSERT_EXCEPTION(false, DAVS,"fail");
    if (l1_epr_update_kind_ != LINEAR) { // :(
        // fetch dx_R and W_R
        SmartPtr<const SymMatrix> W_R = IpData().W();

        // fetch Wx
        SmartPtr<const CompoundSymMatrix> WC = dynamic_cast<const CompoundSymMatrix *>(GetRawPtr(W_R));
        SmartPtr<const SymMatrix> Wx = dynamic_cast<const SymMatrix *>(GetRawPtr(WC->GetComp(0, 0)));
        SmartPtr<Vector> dxx = dx->MakeNew();  // vessel for the final product
        dxx->Set(0.);
        // fetch ds
        SmartPtr<Vector> dss = ds->MakeNew();  // vessel for the final product
        dss->Set(0.);
        // fetch sigma_x
        SmartPtr<const Vector> sigma_R = IpCq().curr_sigma_x();

        SmartPtr<const CompoundVector> sigmaC = dynamic_cast<const CompoundVector *>(GetRawPtr(sigma_R));
        SmartPtr<const Vector> sigma_x = sigmaC->GetComp(0);
        // fetch sigma_s
        SmartPtr<const Vector> sigma_s = IpCq().curr_sigma_s();
        //


        if (l1_epr_update_kind_ == QUAD) {
            Tmp_Sigma_x().SetDiag(*sigma_x);
        }
        else
        {
            SmartPtr<Vector> sigma_dummy_0 = sigma_x->MakeNew();
            sigma_dummy_0->Set(0.);
            Tmp_Sigma_x().SetDiag(*sigma_dummy_0);
        }


        Tmp_Wx_Sigma_x().SetTerm(0, 1.0, *Wx);
        Tmp_Wx_Sigma_x().SetTerm(1, 1.0, Tmp_Sigma_x());
        Tmp_Wx_Sigma_x().MultVector(0.5, *dx, 0.0, *dxx);
        dxWxdx = dxx->Dot(*dx);
        dss->Copy(*ds);
        dss->Scal(0.5);
        dss->ElementWiseMultiply(*sigma_s);
        // put condition here
        if (l1_epr_update_kind_ == QUADNOSIGNA)
        {
            dss->Set(0.0);
        }
        dsWsdx = dss->Dot(*ds);
        // sigma_w
        if ((dxWxdx + dsWsdx) > 0.)
        {
            sigW = 1;
        }
    }
    DBG_PRINT((0, "denominator %8.2e\n", (1.0-l1_epr_epsi_) * ck - ppnp));

    // grad_barrier parts
    Number dphidx = grad_phi_x->Dot(*dx);
    Number dphids = grad_phi_s->Dot(*ds);
    DBG_PRINT((0, "directional %8.2e\n", dphidx + dphids));
    rho_trial = ((dxWxdx + dsWsdx) * sigW + dphidx + dphids)/((1.0-l1_epr_epsi_) * ck - ppnp);
    l1_epr_suff_feasib_update_ = ((1.0 - l1_epr_epsi_) * ck - ppnp) > 0;
    DBG_PRINT((2,"(1-%e)||c|| - (p+n)e=%e \n", l1_epr_epsi_,(1.0-l1_epr_epsi_) * ck - ppnp));
    trial_rho_cache_.AddCachedResult(rho_trial, tdeps, sdeps);  // ToDo have this working properly
    return rho_trial;

}

SumSymMatrix &L1ExactPenaltyRhoUpdater::Tmp_Wx_Sigma_x() {
    if (!IsValid(Hx_Sigma_x_))
    {
        // Create the new Matrix
        Hx_Sigma_x_ = Hx_p_Sigma_x_space_->MakeNewSumSymMatrix();
    }
    return *Hx_Sigma_x_;
}

void L1ExactPenaltyRhoUpdater::SetMatrixSpaces()
{
    // Make new matrix by using the original h_space
    auto* IpL1NLP = dynamic_cast<L1ExactPenaltyRestoIpoptNLP*>(&IpNLP());
    SmartPtr<const VectorSpace> orig_x_space;
    SmartPtr<const VectorSpace> orig_c_space;
    SmartPtr<const VectorSpace> orig_d_space;
    SmartPtr<const VectorSpace> orig_x_l_space;
    SmartPtr<const MatrixSpace> orig_px_l_space;
    SmartPtr<const VectorSpace> orig_x_u_space;
    SmartPtr<const MatrixSpace> orig_px_u_space;
    SmartPtr<const VectorSpace> orig_d_l_space;
    SmartPtr<const MatrixSpace> orig_pd_l_space;
    SmartPtr<const VectorSpace> orig_d_u_space;
    SmartPtr<const MatrixSpace> orig_pd_u_space;
    SmartPtr<const MatrixSpace> orig_jac_c_space;
    SmartPtr<const MatrixSpace> orig_jac_d_space;


    IpL1NLP->OrigIpNLP().GetSpaces(orig_x_space,
                                   orig_c_space,
                                   orig_d_space,
                                   orig_x_l_space,
                                   orig_px_l_space,
                                   orig_x_u_space,
                                   orig_px_u_space,
                                   orig_d_l_space,
                                   orig_pd_l_space,
                                   orig_d_u_space,
                                   orig_pd_u_space,
                                   orig_jac_c_space,
                                   orig_jac_d_space,
                                   original_h_space_);
    // Create Diagonal space for sigma_x
    Sigma_x_space_ = new DiagMatrixSpace(original_h_space_->NRows());
    // Create SumSymMaxtrix space for W + sigma_x
    Hx_p_Sigma_x_space_ = new SumSymMatrixSpace(original_h_space_->NRows(), 2);

}

DiagMatrix &L1ExactPenaltyRhoUpdater::Tmp_Sigma_x()
{
    if (!IsValid(Sigma_x_))
    {
        Sigma_x_ = Sigma_x_space_->MakeNewDiagMatrix();
    }
    return *Sigma_x_;
}

bool L1ExactPenaltyRhoUpdater::UpdateRhoTrial() {
    Number trial_rho = ComputeRhoTrial();
    Number old_rho = L1EPRAddData().GetCurrentRho();
    Number new_rho = old_rho;

    if (l1_epr_update_kind_ == CONST)
    {
        L1EPRAddData().SetRhoTrial(new_rho);
        return false;
    } // if fixed

    l1_epr_has_changed_ = true;

    if (old_rho < trial_rho)
    {
        new_rho = trial_rho + 1.;
        new_rho = new_rho < l1_epr_max_rho * 10. ? new_rho : l1_epr_max_rho * 10.; // Put a maximum value on this x10 times.


        IpData().Append_info_string(" rhoT");
    }
    else if (!l1_epr_suff_feasib_update_ && old_rho < l1_epr_max_rho)
    {
        new_rho = l1_epr_max_rho < old_rho * 2. ? l1_epr_max_rho : old_rho * 2.;
        IpData().Append_info_string(" rhoF");
    }
    else
    {
        l1_epr_has_changed_ = false;
    }
    L1EPRAddData().SetRhoTrial(new_rho);
    L1EPRAddData().SetRhoStatus(l1_epr_has_changed_);
    Jnlst().Printf(J_DETAILED, J_MAIN, "Rho val %7.2e\n", new_rho);
    return l1_epr_has_changed_;
    }

void L1ExactPenaltyRhoUpdater::UpdateRhoAction()
{
    L1EPRAddData().AcceptRhoTrial();
    // Make this optional.
    // If we use a filter AND if the penalty has changed
    if( resto_lsacceptor_option_ == "filter" && l1_epr_has_changed_){
        auto ip_flsa_ = dynamic_cast<FilterLSAcceptor*>(GetRawPtr(ip_bls_acceptor_));
        ip_flsa_->Reset();
    }
// reset filter!;
}

L1ExactPenaltyRestoData &L1ExactPenaltyRhoUpdater::L1EPRAddData()
{
    auto retval = dynamic_cast<L1ExactPenaltyRestoData*>(&IpData().AdditionalData());
    return *retval;
}




}

