//
// Created by David on 1/9/2021.
//

#include "IpL1ExactPenaltyRho.hpp"
#include "IpCompoundVector.hpp"

Ipopt::Number Ipopt::L1ExactPenaltyRho::ComputeRhoTrial()
{

    if (l1exactpenalty_rho_type_ == CONST){return l1exactpenalty_rho0_;} // if fixed

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
    curr_rho_cache_.GetCachedResult(rho_trial, tdeps, sdeps);

    SmartPtr<const CompoundVector> xC = dynamic_cast<const CompoundVector*>(GetRawPtr(x_resto));

    // fetch grad_f_barrier
    SmartPtr<const Vector> grad_barrier_xR = curr_grad_barrier_obj_x();
    SmartPtr<const CompoundVector> grad_barrier_xRC = dynamic_cast<const CompoundVector*>(GetRawPtr(grad_barrier_xR));
    SmartPtr<const Vector> grad_phi_x = grad_barrier_xRC->GetComp(0);
    SmartPtr<const Vector> grad_phi_s = curr_grad_barrier_obj_s();
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
    if (!(l1exactpenalty_rho_type_ == LINEAR)) { // :(
        // fetch dx_R and W_R
        SmartPtr<const SymMatrix> W_R = ip_data_->W();

        // fetch Wx
        SmartPtr<const CompoundSymMatrix> WC = dynamic_cast<const CompoundSymMatrix *>(GetRawPtr(W_R));
        SmartPtr<const SymMatrix> Wx = dynamic_cast<const SymMatrix *>(GetRawPtr(WC->GetComp(0, 0)));
        SmartPtr<Vector> dxx = dx->MakeNew();  // vessel for the final product
        dxx->Set(0.);
        // fetch ds
        SmartPtr<Vector> dss = ds->MakeNew();  // vessel for the final product
        dss->Set(0.);
        // fetch sigma_x
        SmartPtr<const Vector> sigma_R = curr_sigma_x();

        SmartPtr<const CompoundVector> sigmaC = dynamic_cast<const CompoundVector *>(GetRawPtr(sigma_R));
        SmartPtr<const Vector> sigma_x = sigmaC->GetComp(0);
        // fetch sigma_s
        SmartPtr<const Vector> sigma_s = curr_sigma_s();
        //
        SmartPtr<DiagMatrixSpace> diag_space_x = new DiagMatrixSpace(Wx->NRows());
        SmartPtr<SumSymMatrixSpace> ss_space = new SumSymMatrixSpace(Wx->NRows(), 2);
        SmartPtr<const SymMatrixSpace> Wx_space = Wx->OwnerSymMatrixSpace();
        ss_space->SetTermSpace(0, *Wx_space);
        ss_space->SetTermSpace(1, *diag_space_x);

        SmartPtr<DiagMatrix> diag_sigma_x = diag_space_x->MakeNewDiagMatrix();

        SmartPtr<Vector> sigma_dummy = sigma_x->MakeNewCopy();
        sigma_dummy->Set(0.0);
        if (l1exactpenalty_rho_type_ == QUAD) {
            diag_sigma_x->SetDiag(*sigma_x);
        }
        diag_sigma_x->SetDiag(*sigma_dummy);
        SmartPtr<SumSymMatrix> Wx_Sigma = ss_space->MakeNewSumSymMatrix();
        Wx_Sigma->SetTerm(0, 1.0, *Wx);
        Wx_Sigma->SetTerm(1, 1.0, *diag_sigma_x);
        Wx_Sigma->MultVector(0.5, *dx, 0.0, *dxx);
        dxWxdx = dxx->Dot(*dx);
        dss->Copy(*ds);
        dss->Scal(0.5);
        dss->ElementWiseMultiply(*sigma_s);
        // put condition here
        if (l1exactpenalty_rho_type_ == QUADNOSIGMA) { // :(
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
    curr_rho_cache_.AddCachedResult(rho_trial, tdeps, sdeps);  // ToDo have this working properly
    return rho_trial;

    return 0;
}
