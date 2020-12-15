//
// Created by David on 10/1/2020.
//

#include "IpL1ExactPenaltyRestoIpoptNlp.hpp"
#include "IpSumSymMatrix.hpp"

namespace Ipopt
{
    IpL1ExactPenaltyRestoIpoptNLP::IpL1ExactPenaltyRestoIpoptNLP(
            IpoptNLP &orig_ip_nlp, IpoptData &orig_ip_data,
            IpoptCalculatedQuantities &orig_ip_cq)
            : RestoIpoptNLP(orig_ip_nlp, orig_ip_data, orig_ip_cq)
    {
        SmartPtr<OrigIpoptNLP> mynlp = static_cast<OrigIpoptNLP*>(&OrigIpNLP());
        SmartPtr<IpoptNLP> mynlp2 = &OrigIpNLP();
        Number mycpu = mynlp->TotalFunctionEvaluationCpuTime();
        std::cout << mycpu << std::endl;
        RestoIpoptNLP::f_evals();
    }

    IpL1ExactPenaltyRestoIpoptNLP::~IpL1ExactPenaltyRestoIpoptNLP() noexcept
    { }

    Number IpL1ExactPenaltyRestoIpoptNLP::f(
            const Vector& x,
            Number        rho
            )
    {
        DBG_START_METH("IpL1ExactPenaltyRestoIpoptNLP::f",
                       dbg_verbosity);
        // fact_f * (pcTe + ncTe + pdT*e + ndT*e) + fact_c * f(x)
        Number ret = 0.0;
        Number fact_f = 1.0;
        Number fact_c = 1.0;
        // Get the factors
        if (l1_exact_penalty_type_ == OBJECTIVE_INV)
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

    SmartPtr<const Vector> IpL1ExactPenaltyRestoIpoptNLP::grad_f(
            const Vector &x,
            Number rho
    )
    {
        Number fact_f = 1.0;
        Number fact_c = 1.0;
        if (l1_exact_penalty_type_ == OBJECTIVE_INV)
        {
            fact_f = 1/rho;
        }
        else
        {
            fact_c = rho;
        }

        SmartPtr<Vector> retPtr = x.MakeNew();
        retPtr->Set(fact_c);

        const CompoundVector* c_vec_in = static_cast<const CompoundVector*>(&x);
        SmartPtr<const Vector> x_only_in = c_vec_in->GetComp(0);

        SmartPtr<const Vector> orig_gf = OrigIpNLP().grad_f(*x_only_in);

        // Get x component from retptr, copy orig_gf, scale
        //       CompoundVector* c_vec = static_cast<CompoundVector*>(GetRawPtr(retPtr));
        CompoundVector* c_vec = static_cast<CompoundVector*>(GetRawPtr(retPtr));
        SmartPtr<Vector> x_only = c_vec->GetCompNonConst(0);
        x_only->Copy(*orig_gf);
        x_only->Scal(fact_f);
        return ConstPtr(retPtr);
    }

    SmartPtr<const SymMatrix> IpL1ExactPenaltyRestoIpoptNLP::h(const Vector &x,
                                                               Number obj_factor,
                                                               const Vector &yc,
                                                               const Vector &yd,
                                                               Number rho)
    {
        Number fact_f = 1.0;
        Number fact_c = 1.0;
        if (l1_exact_penalty_type_ == OBJECTIVE_INV)
        {
            fact_f = 1/rho;
        }
        else
        {
            fact_c = rho;
        }

        // Get the x_only part
        const CompoundVector* c_in = static_cast<const CompoundVector*>(&x);
        SmartPtr<const CompoundVector> x_in = c_in->GetComp(0);

        // yc and yd (get the 0th)
        const CompoundVector* yc_c = static_cast<const CompoundVector*>(&yc);
        SmartPtr<const CompoundVector> yc_0 = yc_c->GetComp(0);
        const CompoundVector* yd_c = static_cast<const CompoundVector*>(&yd);
        SmartPtr<const CompoundVector> yd_0 = yd_c->GetComp(0);

        // calculate the original hessian

        SmartPtr<const SymMatrix> h_orig = OrigIpNLP().h(*x_in, fact_f, *yc_0, *yd_0);

        // set the compound matrix from the original hessian

        const CompoundSymMatrixSpace* h_comp_space =
                static_cast<const CompoundSymMatrixSpace*>(GetRawPtr(RestoIpoptNLP::HessianMatrixSpace()));

        SmartPtr<CompoundSymMatrix> retPtr = h_comp_space->MakeNewCompoundSymMatrix();
        SmartPtr<Matrix> h_sum_mat = retPtr->GetCompNonConst(0, 0);
        SumSymMatrix* h_sum = static_cast<SumSymMatrix*>(GetRawPtr(h_sum_mat));
        h_sum->SetTerm(0, 1.0, *h_orig);

        return GetRawPtr(retPtr);

    }

    SmartPtr<const SymMatrix> IpL1ExactPenaltyRestoIpoptNLP::uninitialized_h()
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

}
