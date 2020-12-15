//
// Created by David on 9/23/2020.
//

#ifndef IPOPT_NEW_IPL1EXACTPENALTYRESTO_HPP
#define IPOPT_NEW_IPL1EXACTPENALTYRESTO_HPP

#include "IpRestoPhase.hpp"
#include "IpIpoptAlg.hpp"
#include "IpEqMultCalculator.hpp"

namespace Ipopt
{
/** Restoration Phase that minimizes the penalty of the  l-1 norm of the
 * constraint violation and the objective function - using the interior point
 * method (Ipopt)
 */
class L1ExactPenaltyRestorationPhase: public RestorationPhase
{
public:
    /**@name Constructor/Destructor */
    //@{
    /** Constructor, takes algorithm objects.
     *
     * The resto_alg object is the restoration phase Ipopt algorithm object.
     * The eq_mult_calculator is used to reinitialize the equality constraint
     * multipliers after restoration phase is finished, unless it is NULL.
     */
     L1ExactPenaltyRestorationPhase(
             IpoptAlgorithm& resto_alg,
             const SmartPtr<EqMultiplierCalculator>& eq_mult_calculator
             );

     /** Destructor */
     virtual ~L1ExactPenaltyRestorationPhase();
     //@}

     virtual bool InitializeImpl(
             const OptionsList& options,
             const std::string& prefix
             );

     static void RegisterOptions(
             SmartPtr<RegisteredOptions> roptions
             );

protected:
    virtual bool PerformRestoration();

private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     *
     * These methods are not implemented and we do not want the compiler to
     * implement them for us, so we declare them private and do not define them.
     * This ensures that they will not be implicitly called/created.
     */
     //@{
     /** Default Constructor */
     L1ExactPenaltyRestorationPhase();

     /** Copy Constructor */
     L1ExactPenaltyRestorationPhase(
             const L1ExactPenaltyRestorationPhase&
             );

     /** Default Assignment Operator */
     void operator=(
             const L1ExactPenaltyRestorationPhase&
             );
     //@}

     /**@name Strategy Objects */
     //@{
     SmartPtr<IpoptAlgorithm> resto_alg_;
     SmartPtr<EqMultiplierCalculator> eq_mult_calculator_;
     //@}

     /** Copy of original options, which is required to initialize the
      * Ipopt algorithm strategy before the restoration phase is started.
      */
      SmartPtr<OptionsList> resto_options_;

      /** @name Algorithmic parameters */
      //@{
      Number constr_mult_reset_threshold_;

      /** Maximal allowed value of a bound multiplier after restoration phase.
       */
       Number bound_mult_reset_threshold_;

       /** Indicates whether the problem can be expected to be infeasible.
        *
        * This will request the object to set the kappa_resto to a small value
        * for the first time rthe restoration phase is called.
        */
        bool expect_infeasible_problem_;

        /** Constraint violation tolerance */
        Number constr_viol_tol;

        /** Pprimal infeasibility tolerance for declaring failure of restoration
         * phase when the non-regular termination test are met.
         */
         Number resto_failure_feasibility_threshold_;
         //@}

         /** Counter for the number of times that the PerformRestoration is
          * called.
          */
          Index count_restorations_;

          /** @name Auxiliary methods */
          //@{
          /** Method for computing "primal-dual" step in bound multipliers,
           * given step in slacks;
           */
           void ComputeBoundMultiplierStep(
                   Vector& delta_z,
                   const Vector& curr_z,
                   const Vector& curr_slack,
                   const Vector& trial_slack
                   );
          //@}

};
} // namespace Ipopt

#endif //IPOPT_NEW_IPL1EXACTPENALTYRESTO_HPP
