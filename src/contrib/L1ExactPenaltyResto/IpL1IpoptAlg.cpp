//
// Created by David on 1/6/2021.
//

#include "IpL1IpoptAlg.hpp"
#include "IpoptConfig.h"
#include "IpJournalist.hpp"
#include "IpRestoPhase.hpp"
#include "IpOrigIpoptNLP.hpp"



namespace Ipopt
{


IpL1IpoptAlg::IpL1IpoptAlg(
        const SmartPtr<SearchDirectionCalculator> &search_dir_calculator,
        const SmartPtr<LineSearch> &line_search,
        const SmartPtr<MuUpdate> &mu_update,
        const SmartPtr<ConvergenceCheck> &conv_check,
        const SmartPtr<IterateInitializer> &iterate_initializer,
        const SmartPtr<IterationOutput> &iter_output,
        const SmartPtr<HessianUpdater> &hessian_updater,
        const SmartPtr<EqMultiplierCalculator> &eq_multiplier_calculator)
        : search_dir_calculator_(search_dir_calculator),
        line_search_(line_search),
        mu_update_(mu_update),
        conv_check_(conv_check),
        iterate_initializer_(iterate_initializer),
        iter_ouput_(iterate_initializer),
        hessian_updater_(hessian_updater),
        eq_multiplier_calculator_(eq_multiplier_calculator)
{
    DBG_START_METH("IpL1IpoptAlg::IpL1IpoptAlg", dbg_verbosity);
    DBG_ASSERT(IsValid(search_dir_calculator_));
    DBG_ASSERT(IsValid(line_search_));
    DBG_ASSERT(IsValid(mu_update_));
    DBG_ASSERT(IsValid(conv_check_));
    DBG_ASSERT(IsValid(iterate_initializer_));
    DBG_ASSERT(IsValid(iter_output_));
    DBG_ASSERT(IsValid(hessian_updater_));
}

void IpL1IpoptAlg::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
{}

static bool copyright_message_printed = true;

bool IpL1IpoptAlg::InitializeImpl(const OptionsList &options,
                                  const std::string &prefix)
{

    DBG_START_METH("IpL1IpoptAlg::RegisterOptions", dbg_verbosity);
    SmartPtr<const OptionsList> my_options;
    my_options = &options;
    copyright_message_printed = true;

    options.GetStringValue("linear_solver", linear_solver_, prefix);

    bool retvalue = IpData().Initialize(Jnlst(), *my_options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the IpIpoptData object failed to initialize.");
    retvalue = IpCq().Initialize(Jnlst(), *my_options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "The IpIpoptCalculatedQuantities object failed to initialize.");
    retvalue = IpNLP().Initialize(Jnlst(), *my_options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "The IpIpoptNLP object failed to initialize.");

    retvalue = iterate_initializer_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "The iterate_initializer strategy failed to initialize.");

    retvalue = mu_update_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "The mu_update strategy failed to initialize.");

    retvalue = search_dir_calculator_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "The search_direction_calculator strategy failed to initialize.");

    retvalue = line_search_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "The line_search strategy failed to initialize.");

    retvalue = conv_check_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "The conv_check strategy failed to initialize.");

    retvalue = iter_ouput_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "THE iter_output strategy failed to initialize.");

    retvalue = hessian_updater_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "The hessian_updater strategy failed to initialize.");

    my_options->GetNumericValue("kappa_sigma", kappa_sigma_, prefix);
    if(!my_options->GetBoolValue("recalc_y", recalc_y_, prefix))
    {
        Index enum_init;
        if(my_options->GetEnumValue("hessian_approximation", enum_init, prefix))
        {
            HessianApproximationType hessian_approximation = HessianApproximationType(enum_init);
            if(hessian_approximation == LIMITED_MEMORY)
            {
                recalc_y_ = true;
            }
        }
    }

    if(recalc_y_)
    {
        my_options->GetNumericValue("recalc_y_feas_tol", recalc_y_feas_tol_, prefix);
    }

    if (prefix == "resto.")
    {
        skip_print_problem_stats_ = true;
    }
    else
    {
        skip_print_problem_stats_ = false;
    }

    return  true;
}

SolverReturn IpL1IpoptAlg::Optimize(bool isResto)
{
    DBG_START_METH("IpL1IpoptAlg::Optimize", dbg_verbosity);

    IpData().TimingStats().OverallAlgorithm().Start();
    IpData().ResetCpuStartTime();

    SolverReturn retval = UNASSIGNED;

    // Main section of the algorithm
    try
    {
        IpData().TimingStats().PrintProblemStatistics().Start();
        InitializeIterates();
        IpData().TimingStats().InitializeIterates().End();

        IpData().TimingStats().CheckConvergence().Start();
        ConvergenceCheck::ConvergenceStatus conv_status = conv_check_->CheckConvergence();
        IpData().TimingStats().CheckConvergence().End();

        // Main Loop
        while(conv_status == ConvergenceCheck::CONTINUE)
        {
            IpData().TimingStats().UpdateHessian().Start();
            UpdateHessian();
            IpData().TimingStats().UpdateHessian().End();

            IpData().TimingStats().OutputIteration().Start();
            OutputIteration();
            IpData().ResetInfo();
            IpData().TimingStats().OutputIteration().End();

            bool emergency_mode = false;

            IpData().TimingStats().UpdateBarrierParameter().Start();
            emergency_mode = !UpdateBarrierParameter();
            IpData().TimingStats().UpdateBarrierParameter().End();

            if(!emergency_mode)
            {
                IpData().TimingStats().ComputeSearchDirection().Start();
                emergency_mode = !ComputeSearchDirection();
                IpData().TimingStats().ComputeSearchDirection().End();
            }

            if(emergency_mode)
            {
                bool retval = line_search_->ActivateFallbackMechanism();
                if(retval)
                {
                    Jnlst().Printf(J_WARNING, J_MAIN,
                                   "WARNING: Problem in step computation; switching to emergency mode.\n");
                }
                else
                {
                    Jnlst().Printf(J_ERROR, J_MAIN,
                                   "ERROR: Problem in step computation, but emergency mode cannot be activated.\n");
                    THROW_EXCEPTION(STEP_COMPUTATION_FAILED, "Step computation failed.");
                }
            }

            IpData().TimingStats().ComputeAcceptableTrialPoint().Start();
            ComputeAcceptableTrialPoint();
            IpData().TimingStats().ComputeAcceptableTrialPoint().End();

            IpData().TimingStats().AcceptTrialPoint().Start();
            AcceptTrialPoint();
            IpData().TimingStats().AcceptTrialPoint().End();

            IpData().Set_iter_count(IpData().iter_count() + 1);

            IpData().TimingStats().CheckConvergence().Start();
            conv_status = conv_check_->CheckConvergence();
            IpData().TimingStats().CheckConvergence().End();
        }

        IpData().TimingStats().OutputIteration().Start();
        OutputIteration();
        IpData().TimingStats().OutputIteration().End();

        if(conv_status == ConvergenceCheck::CONVERGED || conv_status == ConvergenceCheck::CONVERGED_TO_ACCEPTABLE_POINT)
        {
            if(IpCq().IsSquareProblem())
            {
                ComputeFeasibilityMultipliers();
            }
        }
        switch (conv_status)
        {
            case ConvergenceCheck::CONVERGED:
                retval = SUCCESS;
                break;
            case ConvergenceCheck::CONVERGED_TO_ACCEPTABLE_POINT:
                retval = STOP_AT_ACCEPTABLE_POINT;
                break;
            case ConvergenceCheck::MAXITER_EXCEEDED:
                retval = MAXITER_EXCEEDED;
                break;
            case ConvergenceCheck::CPUTIME_EXCEEDED:
                retval = CPUTIME_EXCEEDED;
                break;
            case ConvergenceCheck::DIVERGING:
                retval = DIVERGING_ITERATES;
                break;
            case ConvergenceCheck::USER_STOP:
                retval = USER_REQUESTED_STOP;
                break;
            default:
                retval = INTERNAL_ERROR;
                break;

        }
    }
    catch (TINY_STEP_DETECTED& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        retval = STOP_AT_ACCEPTABLE_POINT;
    }
    catch (ACCEPTABLE_POINT_REACHED& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        if(IpCq().IsSquareProblem())
        {
            ComputeFeasibilityMultipliers();
        }
        retval = STOP_AT_ACCEPTABLE_POINT;
    }
    catch (LOCALLY_INFEASIBLE& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        IpData().TimingStats().CheckConvergence().EndIfStarted();
        retval = LOCAL_INFEASIBILITY;
    }
    catch (RESTORATION_CONVERGED_TO_FEASIBLE_POINT& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        retval = FEASIBLE_POINT_FOUND; // Instead of declaring failure.
    }
    catch(RESTORATION_FAILED& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        retval = MAXITER_EXCEEDED;
    }
    catch (RESTORATION_MAXITER_EXCEEDED& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        retval = CPUTIME_EXCEEDED;
    }
    catch (RESTORATION_USER_STOP& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        retval = USER_REQUESTED_STOP;
    }
    catch (STEP_COMPUTATION_FAILED& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        if( IpCq().IsSquareProblem())
        {
            ComputeFeasibilityMultipliers();
        }
        retval = FEASIBLE_POINT_FOUND;
    }
    catch (IpoptNLP::Eval_Error& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        retval = INVALID_NUMBER_DETECTED;
    }
    catch (FEASIBILITY_PROBLEM_SOLVED& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        if(IpCq().IsSquareProblem())
        {
            ComputeFeasibilityMultipliers();
        }
        retval = FEASIBLE_POINT_FOUND;
    }
    catch (TOO_FEW_DOF& exc)
    {
        exc.ReportException(Jnlst(), J_MOREDETAILED);
        IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
        retval = TOO_FEW_DEGREES_OF_FREEDOM;
    }
    catch (INTERNAL_ABORT& exc)
    {
        exc.ReportException(Jnlst());
        retval = INTERNAL_ERROR;
    }

    DBG_ASSERT(retval != UNASSIGNED && "Unknown return code in the algorithm");
    IpData().TimingStats().OverallAlgorithm().End();
    return retval;
}


void IpL1IpoptAlg::UpdateHessian()
{
    hessian_updater_->UpdateHessian();
}

bool IpL1IpoptAlg::UpdateBarrierParameter()
{
    bool retval = mu_update_->UpdateBarrierParameter();

    return retval;
}

bool IpL1IpoptAlg::ComputeSearchDirection()
{
    bool retval = search_dir_calculator_->ComputeSearchDirection();

    return retval;
}

void IpL1IpoptAlg::ComputeAcceptableTrialPoint()
{
    line_search_->FindAcceptableTrialPoint();
}

void IpL1IpoptAlg::OutputIteration()
{
    iter_ouput_->WriteOutput();
}

void IpL1IpoptAlg::InitializeIterates()
{
    bool retval = iterate_initializer_->SetInitialIterates();

}

void IpL1IpoptAlg::AcceptTrialPoint()
{
    if(line_search_->CheckSkippedLineSearch())
    {
        Jnlst().Printf(J_SUMMARY, J_MAIN,
                       "Line search didn't find acceptable trial point.\n");
        return;
    }

    Index adjusted_slack = IpCq().AdjustedTrialSlacks();
    DBG_PRINT((1, "adjusted_slack = %d\n", adjusted_slack));
}
}


