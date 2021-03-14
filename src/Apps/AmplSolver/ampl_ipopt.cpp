// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "AmplTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "IpoptConfig.h"
#include "IpUtils.hpp"

#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

#include <cstring>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <chrono>
#include <ratio>


int main(
   int argc,
   char** args)
{
   using namespace Ipopt;

   SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
   app->RethrowNonIpoptException(false);

   // Check if executable is run only to print out options documentation
   if( argc == 2 )
   {
      bool print_options = false;
      std::string print_options_mode("text");
      if( !strcmp(args[1], "--print-options=latex") )
      {
         print_options = true;
         print_options_mode = "latex";
      }
      else if( !strcmp(args[1], "--print-options=doxygen") )
      {
         print_options = true;
         print_options_mode = "doxygen";
      }
      else if( !strcmp(args[1], "--print-options") )
      {
         print_options = true;
      }
      else if( !strcmp(args[1], "--print-latex-options") )
      {
         fprintf(stderr, "ampl_ipopt.cpp: Options --print-latex-options has been replaced by --print-options=latex. Please adjust your call.\n");
         exit(-200);
      }
      if( print_options )
      {
         SmartPtr<OptionsList> options = app->Options();
         options->SetStringValue("print_options_documentation", "yes");
         options->SetStringValue("print_options_mode", print_options_mode);
         app->Initialize("");
         return 0;
      }
   }


   // Call Initialize the first time to create a journalist, but ignore
   // any options file
   ApplicationReturnStatus retval;
   retval = app->Initialize("");
   if( retval != Solve_Succeeded )
   {
      printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
      exit(-100);
   }

   // Add the suffix handler for scaling
   SmartPtr<AmplSuffixHandler> suffix_handler = new AmplSuffixHandler();
   suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source,
                                      AmplSuffixHandler::Number_Type);
   suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source,
                                      AmplSuffixHandler::Number_Type);
   suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source,
                                      AmplSuffixHandler::Number_Type);
   // Modified for warm-start from AMPL
   suffix_handler->AddAvailableSuffix("ipopt_zL_out", AmplSuffixHandler::Variable_Source,
                                      AmplSuffixHandler::Number_Type);
   suffix_handler->AddAvailableSuffix("ipopt_zU_out", AmplSuffixHandler::Variable_Source,
                                      AmplSuffixHandler::Number_Type);
   suffix_handler->AddAvailableSuffix("ipopt_zL_in", AmplSuffixHandler::Variable_Source,
                                      AmplSuffixHandler::Number_Type);
   suffix_handler->AddAvailableSuffix("ipopt_zU_in", AmplSuffixHandler::Variable_Source,
                                      AmplSuffixHandler::Number_Type);

   // David's mark
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    auto t2 = t1.time_since_epoch();
    //auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(t1);
    //auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    auto oname = std::to_string(t2.count());
    app->Options()->SetStringValue("output_file", oname);


    SmartPtr<TNLP> ampl_tnlp = new AmplTNLP(ConstPtr(app->Jnlst()), app->Options(), args, suffix_handler);

   // Call Initialize again to process output related options
   retval = app->Initialize();
   if( retval != Solve_Succeeded )
   {
      printf("ampl_ipopt.cpp: Error in second Initialize!!!!\n");
      exit(-101);
   }

   // David's mark
   std::ofstream myfile;
   std::string fname = "l1EPR";
   std::string ending = "_timings.txt";
   myfile.open(fname + ending, std::ios::app);
   myfile << std::endl;
   std::copy(args + 1, args + argc, std::ostream_iterator<char*>(myfile, " ")); // print the name of the problem

   myfile << "\t" << std::fixed << t2.count();
   Index n, m, nnzJ, nnzH;
   Ipopt::TNLP::IndexStyleEnum index_style;
   ampl_tnlp->get_nlp_info(n, m, nnzJ, nnzH, index_style);
   myfile << "\tn\t" << n << "\tm\t" << m;
   myfile.close();

   const int n_loops = 1; // make larger for profiling
   for( Index i = 0; i < n_loops; i++ )
   {
      retval = app->OptimizeTNLP(ampl_tnlp);
   }

   // David's mark
   myfile.open(fname + ending, std::ios::app);
   auto app_data = app->IpoptDataObject();
   auto app_cqs = app->IpoptCQObject();
   auto it_count = app_data->iter_count();
   auto cpu_timing = app_data->TimingStats().OverallAlgorithm().TotalCpuTime();
   std::cout.precision(8);
   myfile << "\tIter\t" << it_count;
   myfile << "\tCPUs\t" << std::scientific << cpu_timing;
   myfile << "\tf:\t" << app_cqs->curr_f();
   myfile << "\tinf_pr:\t" << app_cqs->curr_primal_infeasibility(NORM_MAX);
   myfile << "\tinf_du:\t" << app_cqs->curr_dual_infeasibility(NORM_MAX);
   myfile << "\tSTAT:\t" << retval;
   myfile.close();

   // finalize_solution method in AmplTNLP writes the solution file

   return 0;
}

