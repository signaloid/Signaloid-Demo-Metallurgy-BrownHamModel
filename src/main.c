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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <inttypes.h>
#include <uxhw.h>
#include "utilities.h"
#include "kernel.h"
#include "common.h"

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
	CommandLineArguments        arguments;
	double                      output;
	double *                    monteCarloOutputSamples = NULL;
	clock_t                     start;
	clock_t                     end;
	double                      cpuTimeUsedInSeconds    = 0.0;
	double                      benchmarkOutput         = 0.0;
	double                      inputDistributions[kInputDistributionIndexMax];
	double                      outputVariables[kOutputDistributionIndexMax];
	const char *                applicationDescription = "Precipitate cutting dislocation model from Brown and Ham";
	const char *const           inputVariableNames[kInputDistributionIndexMax]      = { "b", "G", "gamma", "M", "phi", "Rs" };
	const char *const           outputVariableNames[kOutputDistributionIndexMax]    = {
		[kOutputDistributionIndexSigma] = "sigmaCMpa",
	};
	kOutputVariableTypeIndex    outputVariableTypes[kOutputDistributionIndexMax] = {
		[kOutputDistributionIndexSigma] = kOutputVariableTypeDistribution,
	};
	MeanAndVariance             meanAndVariance = { 0 };

	/*
	 *	Get command-line arguments.
	 */
	if (getCommandLineArguments(argc, argv, &arguments) != kCommonConstantReturnTypeSuccess)
	{
		return EXIT_FAILURE;
	}

	/*
	 *	Read input distributions from CSV if input from file is enabled. This
	 *	overwrites `arguments`' `gamma`/`phi`/`Rs`/`G`/`b`/`M` fields with the
	 *	CSV-derived values; input-from-file is not compatible with Monte Carlo
	 *	mode (rejected above, in `getCommandLineArguments()`), so only the
	 *	UxHw kernel ever sees them.
	 */
	if (arguments.common.isInputFromFileEnabled)
	{
		if (readInputDoubleDistributionsFromCSV(
				arguments.common.inputFilePath,
				inputVariableNames,
				inputDistributions,
				kInputDistributionIndexMax
		))
		{
			return EXIT_FAILURE;
		}

		loadInputs(inputDistributions, &arguments);
	}

	/*
	 *	MonteCarlo output samples are used even in the UxHw use case to store
	 *	the result of the single computed sample.
	 */
	monteCarloOutputSamples =
		(double *) checkedMalloc(
			arguments.common.numberOfMonteCarloIterations * sizeof(double),
			__FILE__,
			__LINE__
		);

	/*
	 *	Start timing.
	 */
	if (arguments.common.isTimingEnabled || arguments.common.isBenchmarkingMode)
	{
		start = clock();
	}

	/*
	 *	Dispatch to the mode-specific kernel. The Monte Carlo loop (when
	 *	applicable) lives inside `calculateOutputMonteCarlo`; UxHw mode runs a
	 *	single distributional evaluation inside `calculateOutputUxHw`.
	 */
	bool isSelectedOutputScalar = (arguments.common.outputSelect != kOutputDistributionIndexMax) &&
	                              (outputVariableTypes[arguments.common.outputSelect] == kOutputVariableTypeScalar);

	if (arguments.common.isMonteCarloMode)
	{
		output = calculateOutputMonteCarlo(&arguments, outputVariables, monteCarloOutputSamples);

		/*
		 *	If not doing UxHw version, then approximate the cost of the third phase of
		 *	Monte Carlo (post-processing), by calculating the mean and variance.
		 *	For scalar outputs, the kernel has already written the correct value to
		 *	`outputVariables[outputSelect]`; the sample buffer holds only a single sample
		 *	at index 0, so the mean would be meaningless. This demo has no scalar outputs,
		 *	so this branch always runs in Monte Carlo mode.
		 */
		if (!isSelectedOutputScalar)
		{
			meanAndVariance = calculateMeanAndVarianceOfDoubleSamples(monteCarloOutputSamples, arguments.common.numberOfMonteCarloIterations);
			output          = outputVariables[arguments.common.outputSelect] = meanAndVariance.mean;
		}
	}
	else
	{
		output = calculateOutputUxHw(&arguments, outputVariables, monteCarloOutputSamples);
	}

	/*
	 *	Stop timing.
	 */
	if (arguments.common.isTimingEnabled || arguments.common.isBenchmarkingMode)
	{
		end = clock();
		cpuTimeUsedInSeconds = ((double) (end - start)) / CLOCKS_PER_SEC;
	}

	benchmarkOutput = output;

	/*
	 *	If in benchmarking mode, print timing result in a special format:
	 *		(1) Benchmark output (for calculating Wasserstein distance to reference)
	 *		(2) Time in microseconds
	 */
	if (arguments.common.isBenchmarkingMode)
	{
		printf("%lf %" PRIu64 "\n", benchmarkOutput, (uint64_t) (cpuTimeUsedInSeconds * 1000000));
	}
	else
	{
		/*
		 *	For scalar outputs in Monte Carlo mode, present a copy of the args with MC
		 *	disabled and iterations=1 so the common print routines take their scalar
		 *	code paths instead of computing distribution stats over a one-element buffer.
		 *	This demo has no scalar outputs, so `printArguments` is always identical to
		 *	`arguments.common`.
		 */
		CommonCommandLineArguments printArguments = arguments.common;

		if (arguments.common.isMonteCarloMode && isSelectedOutputScalar)
		{
			printArguments.isMonteCarloMode             = false;
			printArguments.numberOfMonteCarloIterations = 1;
		}

		/*
		 *	Print json outputs if in JSON output mode.
		 */
		if (arguments.common.isOutputJSONMode)
		{
			printJSONFormattedOutput(
				&printArguments,
				monteCarloOutputSamples,
				outputVariables,
				outputVariableNames,
				kOutputDistributionIndexMax,
				applicationDescription
			);
		}
		/*
		 *	Else print human-consumable output.
		 */
		else
		{
			printf("Cutting stress (σc) = %le MPa\n", outputVariables[kOutputDistributionIndexSigma]);
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
	 *	Save Monte Carlo outputs in an output file.
	 */
	if (arguments.common.isMonteCarloMode)
	{
		size_t samplesToSave = isSelectedOutputScalar
		                ? 1
		                : arguments.common.numberOfMonteCarloIterations;

		saveMonteCarloDoubleDataToDataDotOutFile(
			monteCarloOutputSamples,
			(uint64_t) (cpuTimeUsedInSeconds * 1000000),
			samplesToSave
		);
	}
	/*
	 *	Save outputs to file if not in Monte Carlo mode and write to file is enabled.
	 */
	else
	{
		if (arguments.common.isWriteToFileEnabled)
		{
			if (writeOutputDoubleDistributionsToCSV(
					arguments.common.outputFilePath,
					outputVariables,
					outputVariableNames,
					kOutputDistributionIndexMax
			))
			{
				fprintf(stderr, "Error: Could not write to output CSV file \"%s\".\n", arguments.common.outputFilePath);
				free(monteCarloOutputSamples);

				return EXIT_FAILURE;
			}
		}
	}

	/*
	 *	Free dynamically-allocated memory.
	 */
	free(monteCarloOutputSamples);

	return EXIT_SUCCESS;
}
