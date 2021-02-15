//
// Created by David on 10/1/2020.
//

#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"
#include "IpSumSymMatrix.hpp"
#include "IpL1ExactPenaltyRestoData.hpp"


namespace Ipopt
{
L1ExactPenaltyRestoIpoptNLP::L1ExactPenaltyRestoIpoptNLP(
        IpoptNLP &orig_ip_nlp,
        IpoptData &orig_ip_data,
        IpoptCalculatedQuantities &orig_ip_cq,
        const SmartPtr<IpoptData>& l1_ip_data)
        :
        RestoIpoptNLP(orig_ip_nlp, orig_ip_data, orig_ip_cq),
        l1_ip_data_(l1_ip_data)
{

}

L1ExactPenaltyRestoIpoptNLP::~L1ExactPenaltyRestoIpoptNLP() noexcept
{ }

void L1ExactPenaltyRestoIpoptNLP::RegisterOptions(
        SmartPtr<RegisteredOptions> roptions
)
{
    roptions->AddStringOption2(
        "l1_objective_type",
        "Choose where to put the rho term phi(x,rho) = f + rho ||c||1 (constraint) vs 1/rho f(x) + ||c||1",
        "constraint",
        "constraint", "phi(x,rho) = f(x) + rho ||c(x)||1",
        "objective_inv", "phi(x,rho) = (1/rho) * f(x) + ||c(x)||1",
        "long desc"
        );
}

bool L1ExactPenaltyRestoIpoptNLP::Initialize(
    const Journalist &jnlst,
    const OptionsList& options,
    const std::string& prefix)
{
    Index l1_int;
    options.GetEnumValue("l1_objective_type", l1_int, prefix);
    l1_epr_objective_type_ = IpL1ExactPenaltyObjectiveType(l1_int);
    return Ipopt::RestoIpoptNLP::Initialize(jnlst, options, prefix);
}

bool L1ExactPenaltyRestoIpoptNLP::l1_epr_inv_objective_type() const
{
    return l1_epr_objective_type_ == OBJECTIVE_INV;
}


Number L1ExactPenaltyRestoIpoptNLP::f(
        const Vector& x,
        Number        rho
        )
{
    DBG_START_METH("L1ExactPenaltyRestoIpoptNLP::f",
                   dbg_verbosity);
    // fact_f * (pcTe + ncTe + pdT*e + ndT*e) + fact_c * f(x)
    Number ret = 0.0;
    Number fact_f = 1.0;
    Number fact_c = 1.0;
    // Get the factors
    if (l1_epr_objective_type_ == OBJECTIVE_INV)
    {
        fact_f = 1./rho;
    }
    else
    {
        fact_c = rho;
    }
    // Turn x into a compound
    const CompoundVector* c_vec = static_cast<const CompoundVector*>(&x);
    // Get x from the compound
    SmartPtr<const Vector> x_only = c_vec->GetComp(0);
    // Get the sumation of penalty variables times factor_c.
    ret = x.Sum() - x_only->Sum();
    ret = ret * fact_c;
    ret += OrigIpNLP().f(*x_only) * fact_f;
    return ret;
}

SmartPtr<const Vector> L1ExactPenaltyRestoIpoptNLP::grad_f(
        const Vector &x,
        Number rho
)
{
    Number fact_f = 1.0;
    Number fact_c = 1.0;
    if (l1_epr_objective_type_ == OBJECTIVE_INV)
    {
        fact_f = 1/rho;
    }
    else
    {
        fact_c = rho;
    }

    SmartPtr<Vector> retPtr = x.MakeNew();
    retPtr->Set(fact_c);

    // Evaluate the gradient of f using only the x part of x_in
    auto c_vec_in = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(c_vec_in);
    SmartPtr<const Vector> x_only_in = c_vec_in->GetComp(0);
    SmartPtr<const Vector> orig_gf = OrigIpNLP().grad_f(*x_only_in);

    // Get x component from retptr, copy orig_gf, scale
    //       CompoundVector* c_vec = static_cast<CompoundVector*>(GetRawPtr(retPtr));
    auto c_vec = dynamic_cast<CompoundVector*>(GetRawPtr(retPtr));
    DBG_ASSERT(c_vec);
    SmartPtr<Vector> x_only = c_vec->GetCompNonConst(0);
    x_only->AddOneVector(fact_f, *orig_gf, 0.0);
    //x_only->Copy(*orig_gf);
    //x_only->Scal(fact_f);
    return ConstPtr(retPtr);
}

SmartPtr<const SymMatrix> L1ExactPenaltyRestoIpoptNLP::h(const Vector &x,
                                                         Number obj_factor,
                                                         const Vector &yc,
                                                         const Vector &yd,
                                                         Number rho)
{
    Number fact_f = 1.0;
    Number fact_c = 1.0;
    if (l1_epr_objective_type_ == OBJECTIVE_INV)
    {
        fact_f = 1/rho;
    }
    else
    {
        fact_c = rho;
    }

    // Get the x_only part
    const CompoundVector* c_in = static_cast<const CompoundVector*>(&x);
    SmartPtr<const Vector> x_in = c_in->GetComp(0);

    // yc and yd (get the 0th)
    const CompoundVector* yc_c = static_cast<const CompoundVector*>(&yc);
    SmartPtr<const Vector> yc_0 = yc_c->GetComp(0);
    const CompoundVector* yd_c = static_cast<const CompoundVector*>(&yd);
    SmartPtr<const Vector> yd_0 = yd_c->GetComp(0);

    // calculate the original hessian

    SmartPtr<const SymMatrix> h_orig = OrigIpNLP().h(*x_in, fact_f, *yc_0, *yd_0);

    // set the compound matrix from the original hessian

    const CompoundSymMatrixSpace* h_comp_space =
            static_cast<const CompoundSymMatrixSpace*>(GetRawPtr(RestoIpoptNLP::HessianMatrixSpace()));

    SmartPtr<CompoundSymMatrix> retPtr = h_comp_space->MakeNewCompoundSymMatrix();
    SmartPtr<Matrix> h_sum_mat = retPtr->GetCompNonConst(0, 0);
    SumSymMatrix* h_sum = static_cast<SumSymMatrix*>(GetRawPtr(h_sum_mat));
    h_sum->SetTerm(0, 1.0, *h_orig);
    h_sum->SetTerm(1, 1., *getL1DiagMatDummy());

    return GetRawPtr(retPtr);

}

SmartPtr<const SymMatrix> L1ExactPenaltyRestoIpoptNLP::uninitialized_h()
{
    SmartPtr<CompoundSymMatrix> retPtr;
    const CompoundSymMatrixSpace* h_comp_space =
            static_cast<const CompoundSymMatrixSpace*>(GetRawPtr(RestoIpoptNLP::HessianMatrixSpace()));
    HessianApproximationType hessian_approximation_ = HessianApproximationType(0);
    if( hessian_approximation_ == LIMITED_MEMORY )
    {
        retPtr = h_comp_space->MakeNewCompoundSymMatrix();
    }
    else
    {
        SmartPtr<const SymMatrix> h_con_orig = RestoIpoptNLP::OrigIpNLP().uninitialized_h();
        retPtr = h_comp_space->MakeNewCompoundSymMatrix();
        SmartPtr<Matrix> h_sum_mat = retPtr->GetCompNonConst(0, 0);
        SmartPtr<SumSymMatrix> h_sum = static_cast<SumSymMatrix*>(GetRawPtr(h_sum_mat));
        h_sum->SetTerm(0, 1.0, *h_con_orig);
    //            SmartPtr<SymMatrix> DR_mat =
    //            h_sum->SetTerm(1, 1.0,);
    }
    return GetRawPtr(retPtr);
}

Number L1ExactPenaltyRestoIpoptNLP::Rho() const
{
    L1ExactPenaltyRestoData* l1_epr_data = dynamic_cast<L1ExactPenaltyRestoData*>(&l1_ip_data_->AdditionalData());
    Number rho0 = l1_epr_data->GetCurrentRho();
    if (l1_epr_inv_objective_type() == OBJECTIVE_INV)
    {
        rho0 = 1.; // It has to be 1 for the RestoResto and Iterate init.
    }

    return rho0;
}

    Number L1ExactPenaltyRestoIpoptNLP::f(const Vector &x) {
        return RestoIpoptNLP::f(x);
    }

    SmartPtr<const Vector>
    L1ExactPenaltyRestoIpoptNLP::grad_f(const Vector &x) {
        return RestoIpoptNLP::grad_f(x);
    }

    SmartPtr<const Vector> L1ExactPenaltyRestoIpoptNLP::c(const Vector &x) {
        return RestoIpoptNLP::c(x);
    }

    SmartPtr<const SymMatrix>
    L1ExactPenaltyRestoIpoptNLP::h(const Vector &x, Number obj_factor,
                                   const Vector &yc, const Vector &yd) {
        return RestoIpoptNLP::h(x, obj_factor, yc, yd);
    }

    SmartPtr<const DiagMatrix> L1ExactPenaltyRestoIpoptNLP::getL1DiagMatDummy()
    {
        if(!IsValid(l1_diag_mat_dum_))
        {

            const CompoundSymMatrixSpace* h_comp_space =
                    static_cast<const CompoundSymMatrixSpace*>(GetRawPtr(RestoIpoptNLP::HessianMatrixSpace()));

            SmartPtr<const MatrixSpace> h_orig_space = h_comp_space->GetCompSpace(0, 0);
            SmartPtr<DiagMatrixSpace> diag_mat_space = new DiagMatrixSpace(h_orig_space->NCols());

            l1_diag_vec_dummy_ = DR_x()->MakeNew();
            l1_diag_vec_dummy_->Set(0.);

            l1_diag_mat_dum_ = diag_mat_space->MakeNewDiagMatrix();
            l1_diag_mat_dum_->SetDiag(*l1_diag_vec_dummy_);
        }
    return ConstPtr(l1_diag_mat_dum_);
    }

}
