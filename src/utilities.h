/*
 *	Copyright (c) 2022–2024, Signaloid.
 *
 *	Permission is hereby granted, free of charge, to any person obtaining a copy
 *	of this software and associated documentation files (the "Software"), to deal
 *	in the Software without restriction, including without limitation the rights
 *	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *	copies of the Software, and to permit persons to whom the Software is
 *	furnished to do so, subject to the following conditions:
 *
 *	The above copyright notice and this permission notice shall be included in all
 *	copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *	SOFTWARE.
 */

#pragma once

#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include "common.h"


#define kDemoSpecificConstantGammaUniformMin                            (0.15)
#define kDemoSpecificConstantGammaUniformMax                            (0.25)
#define kDemoSpecificConstantPhiUniformMin                              (0.3)
#define kDemoSpecificConstantPhiUniformMax                              (0.45)
#define kDemoSpecificConstantRsMixtureFirstGaussianMean                 (1E-8)
#define kDemoSpecificConstantRsMixtureFirstGaussianStandardDeviation    (2E-9)
#define kDemoSpecificConstantRsMixtureSecondGaussianMean                (3E-8)
#define kDemoSpecificConstantRsMixtureSecondGaussianStandardDeviation   (2E-9)
#define kDemoSpecificConstantRsMixtureFirstGaussianWeight               (0.5)
#define kDemoSpecificConstantGUniformMin                                (6E10)
#define kDemoSpecificConstantGUniformMax                                (8E10)
#define kDemoSpecificConstantB                                          (2.54E-10)
#define kDemoSpecificConstantMUniformMin                                (1.9)
#define kDemoSpecificConstantMUniformMax                                (4.1)

typedef enum
{
	kInputDistributionIndexB = 0,
	kInputDistributionIndexG,
	kInputDistributionIndexGamma,
	kInputDistributionIndexM,
	kInputDistributionIndexPhi,
	kInputDistributionIndexRs,
	kInputDistributionIndexMax,
} InputDistributionIndex;

typedef enum
{
	kOutputDistributionIndexSigma = 0,
	kOutputDistributionIndexMax,
} OutputDistributionIndex;

typedef struct CommandLineArguments
{
	CommonCommandLineArguments  common;
	double                      gamma;
	double                      phi;
	double                      Rs;
	double                      G;
	double                      b;
	double                      M;
	/*
	 *	Track which of the (possibly distributional) inputs above were
	 *	pinned to an explicit value from the command line, as opposed to
	 *	being left at their default distributions. The Monte Carlo kernel
	 *	uses these to decide, independently for each of its iterations,
	 *	whether to reuse the fixed value or draw a fresh sample.
	 */
	bool    isGammaOverridden;
	bool    isPhiOverridden;
	bool    isRsOverridden;
	bool    isGOverridden;
	bool    isMOverridden;
} CommandLineArguments;

/**
 *	@brief	Print out command line usage.
 */
void
printUsage(void);

/**
 *	@brief	Set default command-line arguments.
 *
 *	@param	arguments	: Pointer to struct that stores command-line arguments.
 *	@return			: `kCommonConstantReturnTypeSuccess` if successful, else `kCommonConstantReturnTypeError`.
 */
CommonConstantReturnType
setDefaultCommandLineArguments(CommandLineArguments * arguments);

/**
 *	@brief	Get command line arguments.
 *
 *	@param	argc		: Argument count from `main()`.
 *	@param	argv		: Argument vector from `main()`.
 *	@param	arguments	: Pointer to struct to store arguments.
 *	@return			: `kCommonConstantReturnTypeSuccess` if successful, else `kCommonConstantReturnTypeError`.
 */
CommonConstantReturnType
getCommandLineArguments(int argc, char *  argv[], CommandLineArguments *  arguments);

/**
 *	@brief	If input-from-file is enabled, overwrite `arguments`' `gamma`, `phi`,
 *		`Rs`, `G`, `b`, and `M` fields with the values read from the input
 *		CSV file. Otherwise, `arguments` is left unmodified: it already
 *		holds either the command-line-supplied values or the UxHw
 *		distributional defaults set by `getCommandLineArguments()`.
 *
 *	@param	inputDistributions	: Array that holds parameter value(s) read from input file.
 *	@param	arguments		: Pointer to struct that stores command-line arguments.
 */
void
loadInputs(
	double *                inputDistributions,
	CommandLineArguments *  arguments);
