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

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <uxhw.h>
#include "utilities.h"
#include "common.h"


/**
 *	@brief	Computes the output of the precipitate dislocation model from Brown and Ham.
 *
 *	@param	gamma	: `gamma` variable.
 *	@param	phi	: `phi` variable.
 *	@param	Rs	: `Rs` variable.
 *	@param	G	: `G` variable.
 *	@param	b	: `b` variable.
 *	@param	M	: `M` variable.
 *	@return		: The output of the precipitate dislocation model from Brown and Ham.
 */
static double
computeBrownHamModelOutput(
	double	gamma,
	double	phi,
	double	Rs,
	double	G,
	double	b,
	double	M)
{
	/*
	 *                    ⎛    _________________    ⎞
	 *       ⎛ M ⋅ γ  ⎞   ⎜   ╱8.0 ⋅ γ ⋅ φ ⋅ Rs     ⎟
	 *  σ  = ⎜─────── ⎟ ⋅ ⎜  ╱ ───────────────── - φ⎟
	 *   c   ⎝2.0 ⋅ b ⎠   ⎝╲╱  π ⋅ G ⋅ pow(b, 2)    ⎠
	 */
	return ((M * gamma) / (2.0 * b))*(sqrt((8.0 * gamma * phi * Rs) / (M_PI * G * pow(b, 2))) - phi) / 1000000;
}

/*
 *	Precipitate "cutting" dislocation model from Brown and Ham
 *
 *	Overview:
 *	Models an equation for a materials precipitate "cutting" dislocation model from Brown and Ham.
 *
 *	Inputs:
 *	The inputs and their distributions are:
 *	-	`gamma`:	Uniform(0.15, 0.25)
 *	-	`phi`:		Uniform(0.30, 0.45)
 *	-	`Rs`:		Equal mixture of Gaussian(1E-8, 2E-9) and Gaussian(3E-8, , 2E-9)
 *	-	`G`:		Uniform(6E10, 8E10)
 *	-	`b`:		2.54E-10 (i.e., constant)
 *	-	`M`:		Uniform(1.9, 4.1)
 *
 *	The parameter `gamma` is the APB energy with units J/m^2, `phi` is the
 *	precipitate volume fraction, `Rs` is mean particle radius on plane with units m,
 *	`G` is the shear modulus with units Pa, `b` is the magnitude of the Burgers
 *	vector with units m, and `M` is the Taylor factor.
 *
 *	Outputs:
 *	The output is the cutting stress, `σc` where
 *
 *                    ⎛    _________________    ⎞
 *       ⎛ M ⋅ γ  ⎞   ⎜   ╱8.0 ⋅ γ ⋅ φ ⋅ Rs     ⎟
 *  σ  = ⎜─────── ⎟ ⋅ ⎜  ╱ ───────────────── - φ⎟
 *   c   ⎝2.0 ⋅ b ⎠   ⎝╲╱  π ⋅ G ⋅ pow(b, 2)    ⎠
 *
 */
int
main(int argc, char *  argv[])
{
	CommandLineArguments	arguments;
	double			gamma;
	double			phi;
	double			Rs;
	double			G;
	double			b;
	double			M;
	double			sigmaCMpa;
	double			inputDistributions[kInputDistributionIndexMax];
	const char * const	outputVariableNames[kOutputDistributionIndexMax] = {"sigmaCMpa"};
	const char * const	inputVariableNames[kInputDistributionIndexMax] = {"b", "G", "gamma", "M", "phi", "Rs"};
	double			outputVariables[kOutputDistributionIndexMax];
	clock_t			start;
	clock_t			end;
	double			cpuTimeUsedInSeconds;
	double			benchmarkOutput;
	double *		monteCarloOutputSamples = NULL;
	MeanAndVariance		monteCarloOutputMeanAndVariance = {0};

	/*
	 *	Get command-line arguments.
	 */
	if (getCommandLineArguments(argc, argv, &arguments) != kCommonConstantReturnTypeSuccess)
	{
		return EXIT_FAILURE;
	}

	/*
	 *	Read input distributions from CSV if input from file is enabled.
	 */
	if (arguments.common.isInputFromFileEnabled)
	{
		if(readInputDoubleDistributionsFromCSV(
			arguments.common.inputFilePath,
			inputVariableNames,
			inputDistributions,
			kInputDistributionIndexMax))
		{
			return EXIT_FAILURE;
		}
	}

	/*
	 *	Allocate for monteCarloOutputSamples if in Monte Carlo mode.
	 */
	if (arguments.common.isMonteCarloMode)
	{
		monteCarloOutputSamples = (double *) checkedMalloc(
								arguments.common.numberOfMonteCarloIterations * sizeof(double),
								__FILE__,
								__LINE__);
	}

	/*
	 *	Start timing.
	 */
	if (arguments.common.isTimingEnabled || arguments.common.isBenchmarkingMode)
	{
		start = clock();
	}

	/*
	 *	Execute process kernel in a loop. The size of loop is 1 unless in Monte Carlo mode.
	 */
	for (size_t i = 0; i < arguments.common.numberOfMonteCarloIterations; ++i)
	{
		/*
		 *	Load inputs.
		 */
		loadInputs(
			&gamma,
			&phi,
			&Rs,
			&G,
			&b,
			&M,
			inputDistributions,
			&arguments);

		/*
		 *	Print inputs if in verbose mode.
		 */
		if (arguments.common.isVerbose)
		{
			printf("Anti-phase boundary energy (γ)\t\t= %le J/m^2\n", gamma);
			printf("Precipitate volume fraction (φ)\t\t= %le\n", phi);
			printf("Mean particle radius on plane (Rs)\t\t= %le m\n", Rs);
			printf("Shear modulus (G)\t\t= %le Pa\n", G);
			printf("Magnitude of the Burger's vector (b)\t\t= %le m\n", b);
			printf("Taylor factor (M)\t\t= %le\n", M);
		}

		/*
		 *	Compute the cutting stress predicted by the Brown-Ham Model.
		 */
		sigmaCMpa = computeBrownHamModelOutput(
				gamma,
				phi,
				Rs,
				G,
				b,
				M);

		/*
		 *	If in Monte Carlo mode, populate monteCarloOutputSamples
		 */
		if (arguments.common.isMonteCarloMode)
		{
			monteCarloOutputSamples[i] = sigmaCMpa;
		}
		/*
		 *	Else, if in benchmarking mode, populate benchmarkOutput.
		 */
		else if (arguments.common.isBenchmarkingMode)
		{
			benchmarkOutput = sigmaCMpa;
		}
	}

	/*
	 *	If not doing Laplace version, then approximate the cost of the third phase of
	 *	Monte Carlo (post-processing), by calculating the mean and variance.
	 */
	if (arguments.common.isMonteCarloMode)
	{
		monteCarloOutputMeanAndVariance = calculateMeanAndVarianceOfDoubleSamples(
								monteCarloOutputSamples,
								arguments.common.numberOfMonteCarloIterations);
		benchmarkOutput = monteCarloOutputMeanAndVariance.mean;
	}

	/*
	 *	Stop timing and evaluate timing result.
	 */
	if (arguments.common.isTimingEnabled || arguments.common.isBenchmarkingMode)
	{
		end = clock();
		cpuTimeUsedInSeconds = ((double) (end - start)) / CLOCKS_PER_SEC;
	}

	/*
	 *	Set outputs.
	 */
	outputVariables[kOutputDistributionIndexSigma] = sigmaCMpa;

	/*
	 *	If in benchmarking mode, print timing result in a special format:
	 *		(1) Benchmark output (for calculating Wasserstein distance to reference)
	 *		(2) Time in microseconds
	 */
	if (arguments.common.isBenchmarkingMode)
	{
		printf("%lf %" PRIu64 "\n", benchmarkOutput, (uint64_t)(cpuTimeUsedInSeconds * 1000000));
	}
	else
	{
		/*
		 *	Print json outputs if in JSON output mode.
		 */
		if (arguments.common.isOutputJSONMode)
		{
			printJSONFormattedOutput(
				sigmaCMpa,
				cpuTimeUsedInSeconds,
				&arguments);
		}
		/*
		 *	Else print human-consumable output.
		 */
		else
		{
			printf("Cutting stress (σc) = %le MPa\n", sigmaCMpa);
		}

		/*
		 *	Print timing if timing is enabled.
		 */
		if (arguments.common.isTimingEnabled)
		{
			printf("CPU time used: %" SignaloidParticleModifier "lf seconds\n", cpuTimeUsedInSeconds);
		}
	}

	/*
	 *	Save Monte Carlo data to "data.out" if in Monte Carlo mode.
	 */
	if (arguments.common.isMonteCarloMode)
	{
		saveMonteCarloDoubleDataToDataDotOutFile(
			monteCarloOutputSamples,
			(uint64_t)(cpuTimeUsedInSeconds * 1000000),
			arguments.common.numberOfMonteCarloIterations);
	}
	/*
	 *	Save outputs to file if not in Monte Carlo mode and write to file is enabled.
	 */
	else
	{
		if (arguments.common.isWriteToFileEnabled)
		{
			if(writeOutputDoubleDistributionsToCSV(
				arguments.common.outputFilePath,
				outputVariables,
				outputVariableNames,
				kOutputDistributionIndexMax))
			{
				fprintf(stderr, "Error: Could not write to output CSV file \"%s\".\n", arguments.common.outputFilePath);

				return EXIT_FAILURE;
			}
		}
	}

	/*
	 *	Free allocations.
	 */
	if (arguments.common.isMonteCarloMode)
	{
		free(monteCarloOutputSamples);
	}

	return EXIT_SUCCESS;
}
