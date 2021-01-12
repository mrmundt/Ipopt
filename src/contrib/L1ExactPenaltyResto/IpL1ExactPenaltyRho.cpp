//
// Created by David on 1/9/2021.
//

#include "IpL1ExactPenaltyRho.hpp"
#include "IpCompoundVector.hpp"
#include "IpCompoundSymMatrix.hpp"
#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"

namespace Ipopt
{
bool L1ExactPenaltyRho::InitializeImpl(const OptionsList &options,
                                              const std::string &prefix) {
    return false;
}



Number L1ExactPenaltyRho::ComputeRhoTrial()
{

    if (l1_epr_update_kind_ == CONST)
    {return l1exactpenalty_rho0_;} // if fixed

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
    c->AddOneVector(-1.0, *nplusc, 1.0);
    c->AddOneVector(1.0, *pplusc, 1.0);
    dms->AddOneVector(-1.0, *nplusd, 1.0);
    dms->AddOneVector(1.0, *pplusd, 1.0);
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
    Number dxWxdx = 0.0;
    Number dsWsdx = 0.0;
    Number sigW = 0.0;
    //THROW_EXCEPTION(DAVS, "fail\n");
    //ASSERT_EXCEPTION(false, DAVS,"fail");
    if (!(l1_epr_update_kind_ == LINEAR)) { // :(
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
        if (l1_epr_update_kind_ == QUADNOSIGMA) { // :(
            dss->Set(0.0);
        }
        dsWsdx = dss->Dot(*ds);
        // sigma_w
        if ((dxWxdx + dsWsdx) > 0.) {
            sigW = 1;
        }
    }
    DBG_PRINT((0, "denominator %8.2e\n", (1.0-epsilon_l1_) * ck - ppnp));

    // grad_barrier parts
    Number dphidx = grad_phi_x->Dot(*dx);
    Number dphids = grad_phi_s->Dot(*ds);
    DBG_PRINT((0, "directional %8.2e\n", dphidx + dphids));
    rho_trial = ((dxWxdx + dsWsdx) * sigW + dphidx + dphids)/((1.0-epsilon_l1_) * ck - ppnp);
    sufficient_feasibility_ = ((1.0 - epsilon_l1_) * ck - ppnp) > 0;
    DBG_PRINT((2,"(1-%e)||c|| - (p+n)e=%e \n", epsilon_l1_,(1.0-epsilon_l1_) * ck - ppnp));
    trial_rho_cache_.AddCachedResult(rho_trial, tdeps, sdeps);  // ToDo have this working properly
    return rho_trial;

    return 0;
}

SumSymMatrix &L1ExactPenaltyRho::Tmp_Wx_Sigma_x() {
    if (!IsValid(Hx_Sigma_x_))
    {
        // Create the new Matrix
        Hx_Sigma_x_ = Hx_p_Sigma_x_space_->MakeNewSumSymMatrix();
    }
    return *Hx_Sigma_x_;
}

void L1ExactPenaltyRho::SetMatrixSpaces()
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

DiagMatrix &L1ExactPenaltyRho::Tmp_Sigma_x()
{
    if (!IsValid(Sigma_x_))
    {
        Sigma_x_ = Sigma_x_space_->MakeNewDiagMatrix();
    }
    return *Sigma_x_;
}

}

