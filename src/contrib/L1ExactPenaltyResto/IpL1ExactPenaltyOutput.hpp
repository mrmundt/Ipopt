//
// Created by dav0 on 2/12/21.
// We need comments here.

#ifndef IPOPT_L1EPR_IPL1EXACTPENALTYOUTPUT_HPP
#define IPOPT_L1EPR_IPL1EXACTPENALTYOUTPUT_HPP

#include "IpIterationOutput.hpp"
#include "IpOrigIterationOutput.hpp"

namespace Ipopt
{
class L1ExactPenaltyOutput: public IterationOutput
{
public:
    explicit L1ExactPenaltyOutput(
            const SmartPtr<OrigIterationOutput>& l1_orig_iteration_output
            );

    ~L1ExactPenaltyOutput() override;

    bool InitializeImpl(
            const OptionsList& options,
            const std::string& prefix
            ) override;
    void WriteOutput() override;

    L1ExactPenaltyOutput() = delete;
    L1ExactPenaltyOutput(
            const L1ExactPenaltyOutput&
    ) = delete;
    L1ExactPenaltyOutput& operator=(
            const L1ExactPenaltyOutput&
            ) = delete;

private:
    SmartPtr<OrigIterationOutput> orig_iteration_output_;

    bool print_info_string_{true};
    InfPrOutput inf_pr_output_;

    int print_frequency_iter_{3};

    Number  print_frequency_time_{1.};
};
}


#endif //IPOPT_L1EPR_IPL1EXACTPENALTYOUTPUT_HPP
