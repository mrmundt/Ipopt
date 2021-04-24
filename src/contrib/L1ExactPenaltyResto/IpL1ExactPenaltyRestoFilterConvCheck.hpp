//
// Created by David Thierry on 1/3/2021.
//
// Based on :
// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//
//           A Waechter: moved most code to IpRestoConvCheck.cpp 2008-06-24


#ifndef SRC_IPL1EXACTPENALTYRESTOFILTERCONVCHECK_HPP
#define SRC_IPL1EXACTPENALTYRESTOFILTERCONVCHECK_HPP

#include "IpOptErrorConvCheck.hpp"
#include "IpBacktrackingLSAcceptor.hpp"

namespace Ipopt
{
/** This is the implementation of the l1-EP restoration convergence check.
 *  For this one, we do not check the filter of the original problem.
 */
class L1ExactPenaltyRestoFilterConvCheck: public OptimalityErrorConvergenceCheck
{
public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    L1ExactPenaltyRestoFilterConvCheck();

    /** Destructor */
    ~L1ExactPenaltyRestoFilterConvCheck() override = default;

    bool InitializeImpl(
            const OptionsList& options,
            const std::string& prefix
            ) override;

    ConvergenceStatus CheckConvergence(
            bool call_intermediate_callback = true
            ) override;


    static void RegisterOptions(
            SmartPtr<RegisteredOptions> roptions
            );

    /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    *
    * These methods are not implemented
    * and we do not want the compiler to implement them for us, so we
    * declare them private and do not define them. This ensures that
    * they will not be implicitly created/called.
    */
    //@{
    /** Copy Constructor */
    L1ExactPenaltyRestoFilterConvCheck(const L1ExactPenaltyRestoFilterConvCheck&) = delete;
    /** Default Assignment Operator */
    L1ExactPenaltyRestoFilterConvCheck& operator = (const L1ExactPenaltyRestoFilterConvCheck&) = delete;
    //@}
private:


    Index maximum_iters_l1_;
    Number orig_constr_viol_tol_;
    bool first_resto_iter_;
    Index succesive_resto_iter_;
};

}
#endif //SRC_IPL1EXACTPENALTYRESTOFILTERCONVCHECK_HPP
