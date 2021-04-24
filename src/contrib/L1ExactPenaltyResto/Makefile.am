# Copyright (C) 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Authors:  Andreas Waechter        IBM       2007-06-01

noinst_LTLIBRARIES = lib_l1.la

lib_l1_la_SOURCES = \
IpL1ExactPenaltyRestoCalculatedQuantities.cpp \
IpL1ExactPenaltyRestoData.cpp \
IpL1ExactPenaltyRestoFilterConvCheck.cpp \
IpL1ExactPenaltyRestoIpoptNlp.cpp \
IpL1ExactPenaltyRestoPdfSpaceSolver.cpp \
IpL1ExactPenaltyRestoRegOp.cpp \
IpL1ExactPenaltyRhoUpdater.cpp \
IpL1ExactPenaltyLeastSquareMults.cpp \
IpL1ExactPenaltyOutput.cpp \
IpL1IpoptAlg.cpp \
IpRestoL1ExactPenalty.cpp


AM_CPPFLAGS = \
	-I$(srcdir)/../../Common \
	-I$(srcdir)/../../LinAlg \
	-I$(srcdir)/../../LinAlg/TMatrices \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/../../Algorithm \
	-I$(srcdir)/../../Algorithm/LinearSolvers \
	$(IPOPTLIB_CFLAGS)
