//
// Created by David on 1/6/2021.
//

#include "IpL1IpoptAlg.hpp"
#include "IpoptConfig.h"
#include "IpJournalist.hpp"
#include "IpRestoPhase.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpL1ExactPenaltyRestoCalculatedQuantities.hpp"



namespace Ipopt
{


L1IpoptAlg::L1IpoptAlg(
        const SmartPtr<SearchDirectionCalculator> &search_dir_calculator,
        const SmartPtr<LineSearch> &line_search,
        const SmartPtr<MuUpdate> &mu_update,
        const SmartPtr<ConvergenceCheck> &conv_check,
        const SmartPtr<IterateInitializer> &iterate_initializer,
        const SmartPtr<IterationOutput> &iter_output,
        const SmartPtr<HessianUpdater> &hessian_updater,
        const SmartPtr<L1ExactPenaltyRhoUpdater> &l1exactpenalty_rho_updater,
        const SmartPtr<EqMultiplierCalculator> &eq_multiplier_calculator)
        : search_dir_calculator_(search_dir_calculator),
        line_search_(line_search),
        mu_update_(mu_update),
        conv_check_(conv_check),
        iterate_initializer_(iterate_initializer),
        iter_ouput_(iter_output),
        hessian_updater_(hessian_updater),
        eq_multiplier_calculator_(eq_multiplier_calculator),
        l1exactpenalty_rho_updater_(l1exactpenalty_rho_updater)

{
    DBG_START_METH("L1IpoptAlg::L1IpoptAlg", dbg_verbosity);
    DBG_ASSERT(IsValid(search_dir_calculator_));
    DBG_ASSERT(IsValid(line_search_));
    DBG_ASSERT(IsValid(mu_update_));
    DBG_ASSERT(IsValid(conv_check_));
    DBG_ASSERT(IsValid(iterate_initializer_));
    DBG_ASSERT(IsValid(iter_output_));
    DBG_ASSERT(IsValid(hessian_updater_));
}

void L1IpoptAlg::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
{}

static bool copyright_message_printed = true;

bool L1IpoptAlg::InitializeImpl(const OptionsList &options,
                                const std::string &prefix)
{

    DBG_START_METH("L1IpoptAlg::RegisterOptions", dbg_verbosity);
    SmartPtr<OptionsList> my_options;
    my_options = new OptionsList(options);
    bool retvalue = false;

    my_options->SetStringValue("start_with_resto", "no", false);
    my_options->SetStringValue("resto.start_with_resto", "no", false);
    copyright_message_printed = true;

    //options.GetStringValue("linear_solver", linear_solver_, prefix);

    retvalue = IpData().Initialize(Jnlst(), *my_options, prefix);
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

    retvalue = search_dir_calculator_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);
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

    retvalue = l1exactpenalty_rho_updater_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(), *my_options, prefix);

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

SolverReturn L1IpoptAlg::Optimize(bool isResto)
{
    DBG_START_METH("L1IpoptAlg::Optimize", dbg_verbosity);

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
                    THROW_EXCEPTION(L1_STEP_COMPUTATION_FAILED, "Step computation failed.");
                }
            }
            ComputeRhoTrial();

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

            UpdateRhoAction();
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
        IpData().TimingStats().UpdateBarrierParameter().EndIfStarted();
        retval = STOP_AT_TINY_STEP;
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
    catch (L1_STEP_COMPUTATION_FAILED& exc)
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


void L1IpoptAlg::UpdateHessian()
{
    hessian_updater_->UpdateHessian();
}

bool L1IpoptAlg::UpdateBarrierParameter()
{
    bool retval = mu_update_->UpdateBarrierParameter();

    return retval;
}

bool L1IpoptAlg::ComputeSearchDirection()
{
    bool retval = search_dir_calculator_->ComputeSearchDirection();

    return retval;
}

void L1IpoptAlg::ComputeAcceptableTrialPoint()
{
    line_search_->FindAcceptableTrialPoint();
}

void L1IpoptAlg::OutputIteration()
{
    iter_ouput_->WriteOutput();
}

void L1IpoptAlg::InitializeIterates()
{
    bool retval = iterate_initializer_->SetInitialIterates();
    ASSERT_EXCEPTION(retval, FAILED_INITIALIZATION, "Error while obtaining initial iterates.");

}

void L1IpoptAlg::AcceptTrialPoint()
{
    if(line_search_->CheckSkippedLineSearch())
    {
        Jnlst().Printf(J_SUMMARY, J_MAIN,
                       "Line search didn't find acceptable trial point.\n");
        return;
    }

    Index adjusted_slacks = IpCq().AdjustedTrialSlacks();
    DBG_PRINT((1, "adjusted_slack = %d\n", adjusted_slack));
    if(adjusted_slacks > 0)
    {
        IpCq().ResetAdjustedTrialSlacks();
        Jnlst().Printf(J_WARNING, J_MAIN,
                       "In iteration %d, %d Slack too small, adjusting variable bound \n",
                       IpData().iter_count(), adjusted_slacks);
        if(Jnlst().ProduceOutput(J_VECTOR, J_MAIN))
        {
            IpNLP().x_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_x_L");
            IpNLP().x_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_x_U");
            IpNLP().d_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_d_L");
            IpNLP().d_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_d_U");
        }

        SmartPtr<Vector> new_x_l = IpNLP().x_L()->MakeNew();
        IpNLP().Px_L()->TransMultVector(1.0, *IpData().trial()->x(), 0.0, *new_x_l);
        new_x_l->Axpy(-1.0, *IpCq().trial_slack_s_L());

        SmartPtr<Vector> new_x_u = IpNLP().x_U()->MakeNew();
        IpNLP().Px_U()->TransMultVector(1.0, *IpData().trial()->x(), 0.0, *new_x_u);
        new_x_u->Axpy(1.0, *IpCq().trial_slack_x_U());

        SmartPtr<Vector> new_d_l = IpNLP().d_L()->MakeNew();
        IpNLP().Pd_L()->TransMultVector(1.0, *IpData().trial()->s(), 0.0, *new_d_l);
        new_d_l->Axpy(-1.0, *IpCq().trial_slack_s_L());

        SmartPtr<Vector> new_d_u = IpNLP().d_U()->MakeNew();
        IpNLP().Pd_L()->TransMultVector(1.0, *IpData().trial()->s(), 0.0, *new_d_u);
        new_d_u->Axpy(1.0, *IpCq().trial_slack_s_U());

        IpNLP().AdjustVariableBounds(*new_x_l, *new_x_u, *new_d_l, *new_d_u);

        if(Jnlst().ProduceOutput(J_VECTOR, J_MAIN))
        {
            IpNLP().x_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_x_L");
            IpNLP().x_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_x_U");
            IpNLP().d_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_d_L");
            IpNLP().d_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_d_U");
        }

    }

    bool corrected = false;
    Number max_correction;
    SmartPtr<const Vector> new_z_L;
    max_correction = correct_bound_multiplier(*IpData().trial()->z_L(),
                                              *IpCq().trial_slack_x_L(),
                                              *IpCq().trial_compl_x_L(),
                                              new_z_L);

    if(max_correction>0.)
    {
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Some value in z_L becomes too large - maximal correction = &8.2e\n",
                       max_correction);
        corrected = true;
    }
    SmartPtr<const Vector> new_z_U;
    max_correction = correct_bound_multiplier(*IpData().trial()->z_U(),
                                              *IpCq().trial_slack_x_U(),
                                              *IpCq().trial_compl_x_U(),
                                              new_z_U);
    if(max_correction>0.)
    {
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Some value in z_U becomes too large - maximal correction = %8.2e\n",
                       max_correction);
        corrected = true;
    }
    SmartPtr<const Vector> new_v_L;
    max_correction = correct_bound_multiplier(*IpData().trial()->v_L(),
                                              *IpCq().trial_slack_s_L(),
                                              *IpCq().trial_compl_s_L(),
                                              new_v_L);
    if(max_correction > 0.)
    {
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Some value in v_L becomes too large - maximal correction = &8.2e\n",
                       max_correction);
        corrected = true;
    }
    SmartPtr<const Vector> new_v_U;
    max_correction = correct_bound_multiplier(*IpData().trial()->v_U(),
                                              *IpCq().trial_slack_s_U(),
                                              *IpCq().trial_compl_s_U(),
                                              new_v_U);
    if(max_correction > 0.)
    {
        Jnlst().Printf(J_DETAILED, J_MAIN,
                       "Some value in v_U becomes too large - maximal correction = &8.2e\n",
                       max_correction);
        corrected = true;
    }
    SmartPtr<IteratesVector> trial = IpData().trial()->MakeNewContainer();
    trial->Set_bound_mult(*new_z_L, *new_z_U, *new_v_L, *new_v_U);
    IpData().set_trial(trial);

    if(corrected)
    {
        IpData().Append_info_string("z");
    }

    IpData().AcceptTrialPoint();

    if(recalc_y_){
        if(IpData().curr()->y_c()->Dim() + IpData().curr()->y_d()->Dim() == 0)
        {
            recalc_y_ = false;
        }
    }

    if(recalc_y_ && IpCq().curr_constraint_violation() < recalc_y_feas_tol_)
    {
        if(Jnlst().ProduceOutput(J_MOREDETAILED, J_MAIN))
        {
            Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                           "dual infeasibility before least square multiplier update = %e\n",
                           IpCq().curr_dual_infeasibility(NORM_MAX));
        }
        IpData().Append_info_string("y ")
        DBG_ASSERT(IsValid(eq_multiplier_calculator_));
        if(IpData().curr()->y_c()->Dim() + IpData().curr()->y_d()->Dim() > 0)
        {
            SmartPtr<Vector> y_c = IpData().curr()->y_c()->MakeNew();
            SmartPtr<Vector> y_d = IpData().curr()->y_d()->MakeNew();
            bool retval = eq_multiplier_calculator_->CalculateMultipliers(*y_c, *y_d);
            if(retval)
            {
                SmartPtr<const IteratesVector> curr = IpData().curr();
                SmartPtr<IteratesVector> iterates = curr->MakeNewContainer();
                iterates->Set_x(*curr->x());
                iterates->Set_s(*curr->s());
                iterates->Set_z_L(*curr->z_L());
                iterates->Set_z_U(*curr->z_U());
                iterates->Set_v_L(*curr->v_L());
                iterates->Set_v_U(*curr->v_U());
                iterates->Set_y_c(*y_c);
                iterates->Set_y_d(*y_d);
                IpData().set_trial(iterates);
                IpData().AcceptTrialPoint();
            }
            else
            {
                Jnlst().Printf(J_DETAILED, J_MAIN,
                               "Recalculation of y multipliers skipped because eq_mult_calc returned false. \n");
            }
        }
    }

}

void L1IpoptAlg::PrintProblemStatistics()
{
    if(!Jnlst().ProduceOutput(J_SUMMARY, J_STATISTICS))
    {
        return;
    }

    Index nx_tot, nx_only_lower, nx_both, nx_only_upper;
    calc_number_of_bounds(*IpData().curr()->x(),
                          *IpNLP().x_L(),
                          *IpNLP().x_U(),
                          *IpNLP().Px_L(),
                          *IpNLP().Px_U(),
                          nx_tot, nx_only_lower, nx_both, nx_only_upper);
    Index ns_tot, ns_only_lower, ns_both, ns_only_upper;
    calc_number_of_bounds(*IpData().curr()->s(),
                          *IpNLP().d_L(),
                          *IpNLP().d_U(),
                          *IpNLP().Pd_L(),
                          *IpNLP().Pd_U(),
                          ns_tot,
                          ns_only_lower,
                          ns_both,
                          ns_only_upper);

    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "Total number of variables............................: %8d\n", nx_tot);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "                     variables with only lower bounds: %8d\n", nx_only_lower);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "                variables with lower and upper bounds: %8d\n", nx_both);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "                     variables with only upper bounds: %8d\n", nx_only_upper);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "Total number of equality constraints.................: %8d\n", IpData().curr()->y_c()->Dim());
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "Total number of inequality constraints...............: %8d\n", ns_tot);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "        inequality constraints with only lower bounds: %8d\n", ns_only_lower);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "   inequality constraints with lower and upper bounds: %8d\n", ns_both);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "        inequality constraints with only upper bounds: %8d\n\n", ns_only_upper);
}

void L1IpoptAlg::ComputeFeasibilityMultipliers()
{
    DBG_START_METH("L1IpoptAlg::ComputeFeasibilityMultipliers",
                   dbg_verbosity);
    DBG_ASSERT(IpCq().IsSquareProblem());
    if(IsNull(eq_multiplier_calculator_))
    {
        Jnlst().Printf(J_WARNING, J_SOLUTION,
                       "This is a square problem, but multipliers cannot be recomputed at solutions, since no eq_mult_calculator object is available in L1IpoptAlgorithm\n");
        return;
    }

    SmartPtr<IteratesVector> iterates = IpData().curr()->MakeNewContainer();
    SmartPtr<Vector> tmp = iterates->z_L()->MakeNew();
    tmp->Set(0.);
    iterates->Set_z_L(*tmp);
    tmp = iterates->z_U()->MakeNew();
    tmp->Set(0.);
    iterates->Set_z_U(*tmp);
    tmp = iterates->v_L()->MakeNew();
    tmp->Set(0.);
    iterates->Set_v_L(*tmp);
    tmp = iterates->v_U()->MakeNew();
    tmp->Set(0.);
    iterates->Set_v_U(*tmp);
    SmartPtr<Vector> y_c = iterates->y_c()->MakeNew();
    SmartPtr<Vector> y_d = iterates->y_d()->MakeNew();
    IpData().set_trial(iterates);
    IpData().AcceptTrialPoint();
    bool retval = eq_multiplier_calculator_->CalculateMultipliers(*y_c, *y_d);
    if(retval)
    {
        iterates = IpData().curr()->MakeNewContainer();
        iterates->Set_y_c(*y_c);
        iterates->Set_y_d(*y_d);
        IpData().set_trial(iterates);
        IpData().AcceptTrialPoint();
    }
    else
    {
        Jnlst().Printf(J_WARNING, J_SOLUTION,
                       "Cannot recompute multipliers for feasibility problem."
                       "Error in eq_mult_calculator ;)\n");
    }
}

void L1IpoptAlg::calc_number_of_bounds(const Vector &x, const Vector &x_L,
                                       const Vector &x_U,
                                       const Matrix &Px_L,
                                       const Matrix &Px_U, Index &n_tot,
                                       Index &n_only_lower, Index &n_both,
                                       Index &n_only_upper)
{
    DBG_START_METH("L1IpoptAlg::calc_number_of_bounds",
                   dbg_verbosity);
    n_tot = x.Dim();
    SmartPtr<Vector> tmpx = x.MakeNew();
    SmartPtr<Vector> tmpxL = x_L.MakeNew();
    SmartPtr<Vector> tmpxU = x_U.MakeNew();

    tmpxL->Set(-1.);
    tmpxU->Set(2.);
    Px_L.MultVector(1., *tmpxL, 0., *tmpx);
    Px_U.MultVector(1., *tmpxU, 1., *tmpx);

    DBG_PRINT_VECTOR(2, "x-indicator", *tmpx);

    SmartPtr<Vector> tmpx0 = x.MakeNew();
    tmpx0->Set(0.);

    SmartPtr<Vector> tmpx2 = x.MakeNew();
    tmpx2->Set(-1.);
    tmpx2->Axpy(1., *tmpx);
    tmpx2->ElementWiseMax(*tmpx0);

    n_only_upper = (Index) tmpx2->Asum();

    tmpx->Axpy(-2., *tmpx2);

    tmpx2->Copy(*tmpx);
    tmpx2->ElementWiseMax(*tmpx0);

    n_both = (Index) tmpx2->Asum();

    tmpx->Axpy(-1., *tmpx2);
    tmpx->ElementWiseMax(*tmpx);

    n_only_lower = (Index) tmpx->Asum();
}

Number L1IpoptAlg::correct_bound_multiplier(const Vector &trial_z,
                                            const Vector &trial_slack,
                                            const Vector &trial_compl,
                                            SmartPtr<const Vector> &new_trial_z)
{
    DBG_START_METH("L1IpoptAlg::correct_bound_multiplier",
                   dbg_verbosity);

    if(kappa_sigma_ < 1. || trial_z.Dim() == 0)
    {
        new_trial_z = &trial_z;
        return 0.;
    }

    Number mu;
    if(IpData().FreeMuMode())
    {
        mu = IpCq().trial_avrg_compl();
        mu = Min(mu, 1e3);
    }
    else
    {
        mu = IpData().curr_mu();
    }
    DBG_PRINT((1, "mu = &8.2e\n", mu));
    DBG_PRINT_VECTOR(2, "trial_z", trial_z);

    if(trial_compl.Amax() <= kappa_sigma_ * mu && trial_compl.Min() >= 1. / kappa_sigma_ * mu)
    {
        new_trial_z = &trial_z;
        return 0.;
    }
    SmartPtr<Vector> one_over_s = trial_z.MakeNew();
    one_over_s->Copy(trial_slack);
    one_over_s->ElementWiseReciprocal();

    SmartPtr<Vector> step_z = trial_z.MakeNew();
    step_z->AddTwoVectors(kappa_sigma_ * mu, *one_over_s, -1., trial_z, 0.);

    DBG_PRINT_VECTOR(2, "step_z", *step_z);

    Number max_correction_up = Max(0., -step_z->Min());
    if(max_correction_up > 0.)
    {
        SmartPtr<Vector> tmp = trial_z.MakeNew();
        tmp->Set(0.);
        step_z->ElementWiseMin(*tmp);
        tmp->AddTwoVectors(1., trial_z, 1., *step_z, 0.);
        new_trial_z = GetRawPtr(tmp);
    }
    else
    {
        new_trial_z = &trial_z;
    }

    step_z->AddTwoVectors(1./kappa_sigma_ * mu, *one_over_s, -1., *new_trial_z, 0.);

    Number max_correction_low = Max(0., step_z->Max());
    if(max_correction_low > 0.)
    {
        SmartPtr<Vector> tmp = trial_z.MakeNew();
        tmp->Set(0.);
        step_z->ElementWiseMax(*tmp);
        tmp->AddTwoVectors(1., *new_trial_z, 1., *step_z, 0.);
        new_trial_z = GetRawPtr(tmp);
    }

    DBG_PRINT_VECTOR(2, "new_trial_z", *new_trial_z);

    return Max(max_correction_up, max_correction_low);
}

void L1IpoptAlg::print_copyright_message(
        const Journalist& jnlst
)
{
    jnlst.Printf(J_INSUPPRESSIBLE, J_MAIN,
                 "\n******************************************************************************\n"
                 "This program contains Ipopt, a library for large-scale nonlinear optimization.\n"
                 " Ipopt is released as open source code under the Eclipse Public License (EPL).\n"
                 "         For more information visit https://github.com/coin-or/Ipopt\n"
                 "******************************************************************************\n\n");
    copyright_message_printed = true;
}

void L1IpoptAlg::ComputeRhoTrial()
{
    l1exactpenalty_rho_updater_->UpdateRhoTrial();
}

void L1IpoptAlg::UpdateRhoAction()
{
    l1exactpenalty_rho_updater_->UpdateRhoAction();
}

    L1IpoptAlg::~L1IpoptAlg() = default;


}


