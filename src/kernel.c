/*
 *	Copyright (c) 2024-2026, Signaloid.
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

#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <uxhw.h>
#include "kernel.h"
#include "brown-and-ham-uxhw.h"
#include "brown-and-ham-monte-carlo.h"
#include "utilities.h"

/*
 *	Print the resolved inputs for one sample if verbose mode is enabled. In
 *	UxHw mode this runs once, for the single resolved input set; in Monte
 *	Carlo mode `brownHamModelMonteCarloSample()` calls this once per
 *	iteration, for that iteration's freshly-drawn inputs, matching the
 *	original per-iteration verbose logging.
 */
static void
printVerboseInputsIfEnabled(
	CommandLineArguments *  arguments,
	double                  gamma,
	double                  phi,
	double                  Rs,
	double                  G,
	double                  b,
	double                  M)
{
	if (arguments->common.isVerbose)
	{
		printf("Anti-phase boundary energy (γ)\t\t= %le J/m^2\n", gamma);
		printf("Precipitate volume fraction (φ)\t\t= %le\n", phi);
		printf("Mean particle radius on plane (Rs)\t\t= %le m\n", Rs);
		printf("Shear modulus (G)\t\t= %le Pa\n", G);
		printf("Magnitude of the Burger's vector (b)\t\t= %le m\n", b);
		printf("Taylor factor (M)\t\t= %le\n", M);
	}

	return;
}

double
computeBrownHamModelOutput(
	double  gamma,
	double  phi,
	double  Rs,
	double  G,
	double  b,
	double  M)
{
	/*
	 *                    ⎛    _________________    ⎞
	 *       ⎛ M ⋅ γ  ⎞   ⎜   ╱8.0 ⋅ γ ⋅ φ ⋅ Rs     ⎟
	 *  σ  = ⎜─────── ⎟ ⋅ ⎜  ╱ ───────────────── - φ⎟
	 *   c   ⎝2.0 ⋅ b ⎠   ⎝╲╱  π ⋅ G ⋅ pow(b, 2)    ⎠
	 */
	return ((M * gamma) / (2.0 * b)) * (sqrt((8.0 * gamma * phi * Rs) / (M_PI * G * pow(b, 2))) - phi) / 1000000;
}

double
brownHamModelMonteCarloSample(CommandLineArguments * arguments)
{
	double  gamma;
	double  phi;
	double  Rs;
	double  G;
	double  M;

	/*
	 *	Inputs pinned by the command line stay fixed across every Monte
	 *	Carlo iteration and free inputs are redrawn from their default
	 *	distributions on every call, mirroring `setDefaultCommandLineArguments()`.
	 */
	gamma = arguments->isGammaOverridden
	        ? arguments->gamma
	        : UxHwDoubleUniformDist(kDemoSpecificConstantGammaUniformMin, kDemoSpecificConstantGammaUniformMax);

	phi = arguments->isPhiOverridden
	        ? arguments->phi
	        : UxHwDoubleUniformDist(kDemoSpecificConstantPhiUniformMin, kDemoSpecificConstantPhiUniformMax);

	Rs = arguments->isRsOverridden
	        ? arguments->Rs
	        : UxHwDoubleMixture(
		UxHwDoubleGaussDist(kDemoSpecificConstantRsMixtureFirstGaussianMean, kDemoSpecificConstantRsMixtureFirstGaussianStandardDeviation),
		UxHwDoubleGaussDist(kDemoSpecificConstantRsMixtureSecondGaussianMean, kDemoSpecificConstantRsMixtureSecondGaussianStandardDeviation),
		kDemoSpecificConstantRsMixtureFirstGaussianWeight
	        );

	G = arguments->isGOverridden
	        ? arguments->G
	        : UxHwDoubleUniformDist(kDemoSpecificConstantGUniformMin, kDemoSpecificConstantGUniformMax);

	M = arguments->isMOverridden
	        ? arguments->M
	        : UxHwDoubleUniformDist(kDemoSpecificConstantMUniformMin, kDemoSpecificConstantMUniformMax);

	printVerboseInputsIfEnabled(arguments, gamma, phi, Rs, G, arguments->b, M);

	return computeBrownHamModelOutput(gamma, phi, Rs, G, arguments->b, M);
}

double
calculateOutputUxHw(
	CommandLineArguments *  arguments,
	double *                outputVariables,
	double *                monteCarloOutputSamples)
{
	double sigma;

	printVerboseInputsIfEnabled(
		arguments,
		arguments->gamma,
		arguments->phi,
		arguments->Rs,
		arguments->G,
		arguments->b,
		arguments->M
	);

	sigma = computeSigmaUxHw(
		arguments->gamma,
		arguments->phi,
		arguments->Rs,
		arguments->G,
		arguments->b,
		arguments->M
	);

	outputVariables[kOutputDistributionIndexSigma] = sigma;
	monteCarloOutputSamples[0] = sigma;

	return sigma;
}

double
calculateOutputMonteCarlo(
	CommandLineArguments *  arguments,
	double *                outputVariables,
	double *                monteCarloOutputSamples)
{
	double sigma;

	sigma = computeSigmaMonteCarlo(
		arguments,
		arguments->common.numberOfMonteCarloIterations,
		monteCarloOutputSamples
	);

	outputVariables[kOutputDistributionIndexSigma] = sigma;

	return sigma;
}
