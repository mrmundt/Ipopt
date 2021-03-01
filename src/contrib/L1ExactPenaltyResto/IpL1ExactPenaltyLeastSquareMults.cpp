//
// Created by dav0 on 2/20/21.
//

#include "IpL1ExactPenaltyLeastSquareMults.hpp"
#include "IpL1ExactPenaltyRestoData.hpp"
#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"

namespace Ipopt {
#if IPOPT_VERBOSITY > 0
    static const Index dbg_verbosity = 0;
#endif

L1ExactPenaltyLeastSquareMults::L1ExactPenaltyLeastSquareMults(
        Ipopt::AugSystemSolver &augSystemSolver)
        : EqMultiplierCalculator(),
          augsyssolver_(&augSystemSolver)

{

}

bool L1ExactPenaltyLeastSquareMults::InitializeImpl(const OptionsList &options, const std::string &prefix) {
    l1EprAddData_ = dynamic_cast<L1ExactPenaltyRestoData*>(&IpData().AdditionalData());
    DBG_ASSERT(l1EprAddData);
    is_rho_inv_obj_ = dynamic_cast<L1ExactPenaltyRestoIpoptNLP*>(&IpNLP())->l1_epr_inv_objective_type();
    return augsyssolver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), options, prefix);
}

bool L1ExactPenaltyLeastSquareMults::CalculateMultipliers(Ipopt::Vector &y_c, Ipopt::Vector &y_d) {
    DBG_START_METH(L1ExactPenaltyLeastSquareMults::CalculateMultipliers,
                   dbg_verbosity);
    Number rho = 1.;
    Number s_factor = 1.0;
    SmartPtr<const SymMatrix> zeroW = IpNLP().uninitialized_h();
    SmartPtr<const Matrix> J_c = IpCq().curr_jac_c();
    SmartPtr<const Matrix> J_d = IpCq().curr_jac_d();
    SmartPtr<const Vector> grad_f = IpCq().curr_grad_f();
    SmartPtr<const Matrix> Px_L = IpNLP().Px_L();
    SmartPtr<const Matrix> Px_U = IpNLP().Px_U();
    SmartPtr<const Matrix> Pd_L = IpNLP().Pd_L();
    SmartPtr<const Matrix> Pd_U = IpNLP().Pd_U();
    SmartPtr<const Vector> z_L = IpData().curr()->z_L();
    SmartPtr<const Vector> z_U = IpData().curr()->z_U();
    SmartPtr<const Vector> v_L = IpData().curr()->v_L();
    SmartPtr<const Vector> v_U = IpData().curr()->v_U();

    // Now, the rhs_x must account for 1/rho.
    SmartPtr<Vector> rhs_x = grad_f->MakeNew();
    rhs_x->Copy(*grad_f);
    SmartPtr<Vector> zDummy;
    if (is_rho_inv_obj_){
        SmartPtr<Vector> zX_only;
        CompoundVector* zD_Compound;
        rho = l1EprAddData_->GetCurrentRho();
        //
        zDummy = z_L->MakeNew();
        zDummy->Copy(*z_L);
        zD_Compound = dynamic_cast<CompoundVector*>(GetRawPtr(zDummy));
        zX_only = zD_Compound->GetCompNonConst(0);
        zX_only->Scal(1/rho);
        Px_L->MultVector(1., *zDummy, -1., *rhs_x);
        //
        zDummy = z_U->MakeNew();
        zDummy->Copy(*z_U);
        zD_Compound = dynamic_cast<CompoundVector*>(GetRawPtr(zDummy));
        zX_only = zD_Compound->GetCompNonConst(0);
        zX_only->Scal(1/rho);
        Px_U->MultVector(-1., *zDummy, 1., *rhs_x);

    } else {
        Px_L->MultVector(1., *z_L, -1., *rhs_x);
        Px_U->MultVector(-1., *z_U, 1., *rhs_x);
    }

    // Scale the d multipliers, if necessary.
    s_factor = is_rho_inv_obj_ ? 1./(l1EprAddData_->GetCurrentRho()) : 1.;

    SmartPtr<Vector> rhs_s = IpData().curr()->s()->MakeNew();
    Pd_L->MultVector(s_factor, *v_L, 0, *rhs_s);
    Pd_U->MultVector(-s_factor, *v_U, 1., *rhs_s);

    SmartPtr<Vector> rhs_c = y_c.MakeNew();
    rhs_c->Set(0.);
    SmartPtr<Vector> rhs_d = y_d.MakeNew();
    rhs_d->Set(0.);

    SmartPtr<Vector> sol_x = rhs_x->MakeNew();
    SmartPtr<Vector> sol_s = rhs_s->MakeNew();

    DBG_PRINT_VECTOR(2, "rhs_x", *rhs_x);
    DBG_PRINT_VECTOR(2, "rhs_s", *rhs_s);
    DBG_PRINT_VECTOR(2, "rhs_c", *rhs_c);
    DBG_PRINT_VECTOR(2, "rhs_d", *rhs_d);

    enum ESymSolverStatus retval;
    Index numberOfEvals = rhs_c->Dim() + rhs_d->Dim();
    bool check_NegEvals = augsyssolver_->ProvidesInertia();
    retval = augsyssolver_->Solve(GetRawPtr(zeroW),
                                  0.,
                                  nullptr,
                                  1.,
                                  nullptr,
                                  1.,
                                  GetRawPtr(J_c),
                                  nullptr,
                                  0.,
                                  GetRawPtr(J_d),
                                  nullptr,
                                  0.,
                                  *rhs_x,
                                  *rhs_s,
                                  *rhs_c,
                                  *rhs_d,
                                  *sol_x,
                                  *sol_s,
                                  y_c,
                                  y_d,
                                  check_NegEvals,
                                  numberOfEvals);

    if(retval != SYMSOLVER_SUCCESS)
    {
        // What is the deal with the original LSM computation failing for SQ problems?
        return false;
    }

    DBG_PRINT_VECTOR(2, "sol_x", *sol_x);
    DBG_PRINT_VECTOR(2, "sol_s", *sol_s);
    DBG_PRINT_VECTOR(2, "sol_c", y_c);
    DBG_PRINT_VECTOR(2, "sol_d", y_d);

    return true;
}

}
