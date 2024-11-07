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
#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <uxhw.h>
#include "utilities.h"
#include "common.h"


void
printUsage(void)
{
	fprintf(stderr, "Example: Precipitate Dislocation Model from Brown and Ham - Signaloid version\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: Valid command-line arguments are:\n");
	fprintf(stderr,
		"\t[-o, --output <Path to output CSV file : str>] (Specify the output file.)\n"
		"\t[-M, --multiple-executions <Number of executions : int> (Default: 1)] (Repeated execute kernel for benchmarking.)\n"
		"\t[-T, --time] (Timing mode: Times and prints the timing of the kernel execution.)\n"
		"\t[-v, --verbose] (Verbose mode: Prints extra information about demo execution.)\n"
		"\t[-b, --benchmarking] (Benchmarking mode: Generate outputs in format for benchmarking.)\n"
		"\t[-j, --json] (Print output in JSON format.)\n"
		"\t[-h, --help] (Display this help message.)\n"
		"\t[-g, --apb-energy <gamma: double> (Default: Uniform(%"SignaloidParticleModifier".2lf, %"SignaloidParticleModifier".2lf))] (Set `gamma` variable.)\n"
		"\t[-p, --precipitate-volume-fraction <phi: double> (Default: Uniform(%"SignaloidParticleModifier".2lf, %"SignaloidParticleModifier".2lf))] (Set `phi` variable.)\n"
		"\t[-R, --mean-particle-radius <Rs: double> (Default: UxHwDoubleMixture(Gauss(%"SignaloidParticleModifier".1le, %"SignaloidParticleModifier".1le), Gauss(%"SignaloidParticleModifier".1le, %"SignaloidParticleModifier".1le), %"SignaloidParticleModifier".1lf))] (Set `Rs` variable.)\n"
		"\t[-G, --shear-modulus <G: double> (Default: Uniform(%"SignaloidParticleModifier".1le, %"SignaloidParticleModifier".1le))] (Set `G` variable.)\n"
		"\t[-B, --burgers-vector <b: double> (Default: %"SignaloidParticleModifier".2le)] (Set `b` variable.)\n"
		"\t[-m, --taylor-factor <M: double> (Default: Uniform(%"SignaloidParticleModifier".1lf, %"SignaloidParticleModifier".1lf))] (Set `M` variable.)\n",
		kDemoSpecificConstantGammaUniformMin,
		kDemoSpecificConstantGammaUniformMax,
		kDemoSpecificConstantPhiUniformMin,
		kDemoSpecificConstantPhiUniformMax,
		kDemoSpecificConstantRsMixtureFirstGaussianMean,
		kDemoSpecificConstantRsMixtureFirstGaussianStandardDeviation,
		kDemoSpecificConstantRsMixtureSecondGaussianMean,
		kDemoSpecificConstantRsMixtureSecondGaussianStandardDeviation,
		kDemoSpecificConstantRsMixtureFirstGaussianWeight,
		kDemoSpecificConstantGUniformMin,
		kDemoSpecificConstantGUniformMax,
		kDemoSpecificConstantB,
		kDemoSpecificConstantMUniformMin,
		kDemoSpecificConstantMUniformMax);
	fprintf(stderr, "\n");

	return;
}

CommonConstantReturnType
setDefaultCommandLineArguments(CommandLineArguments *	arguments)
{
	if (arguments == NULL)
	{
		fprintf(stderr, "Error: The provided pointer to arguments is NULL.\n");

		return kCommonConstantReturnTypeError;
	}

/*
 *	Older GCC versions have a bug which gives a spurious warning for the C universal zero
 *	initializer `{0}`. Any workaround makes the code less portable or prevents the common code
 *	from adding new fields to the `CommonCommandLineArguments` struct. Therefore, we surpress
 *	this warning.
 *
 *	See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53119.
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"

	*arguments = (CommandLineArguments) {
		.common			= (CommonCommandLineArguments) {0},
		.gamma			= UxHwDoubleUniformDist(kDemoSpecificConstantGammaUniformMin, kDemoSpecificConstantGammaUniformMax),
		.phi			= UxHwDoubleUniformDist(kDemoSpecificConstantPhiUniformMin, kDemoSpecificConstantPhiUniformMax),
		.Rs 			= UxHwDoubleMixture(
						UxHwDoubleGaussDist(kDemoSpecificConstantRsMixtureFirstGaussianMean, kDemoSpecificConstantRsMixtureFirstGaussianStandardDeviation),
						UxHwDoubleGaussDist(kDemoSpecificConstantRsMixtureSecondGaussianMean, kDemoSpecificConstantRsMixtureSecondGaussianStandardDeviation),
						kDemoSpecificConstantRsMixtureFirstGaussianWeight),
		.G			= UxHwDoubleUniformDist(kDemoSpecificConstantGUniformMin, kDemoSpecificConstantGUniformMax),
		.b			= kDemoSpecificConstantB,
		.M			= UxHwDoubleUniformDist(kDemoSpecificConstantMUniformMin, kDemoSpecificConstantMUniformMax),
	};

	return kCommonConstantReturnTypeSuccess;
}

CommonConstantReturnType
getCommandLineArguments(int argc, char *  argv[], CommandLineArguments *  arguments)
{
	const char *	gammaArg = NULL;
	const char *	phiArg = NULL;
	const char *	RsArg = NULL;
	const char *	GArg = NULL;
	const char *	bArg = NULL;
	const char *	MArg = NULL;
	const char	kConstantStringUx[] = "Ux";

	if (arguments == NULL)
	{
		fprintf(stderr, "Error: The provided pointer to arguments is NULL.\n");

		return kCommonConstantReturnTypeError;
	}

	if (setDefaultCommandLineArguments(arguments) != kCommonConstantReturnTypeSuccess)
	{
		return kCommonConstantReturnTypeError;
	}

	DemoOption options[] = {
		{ .opt = "g", .optAlternative = "apb-energy", .hasArg = true,.foundArg = &gammaArg,	.foundOpt = NULL },
		{ .opt = "p", .optAlternative = "precipitate-volume-fraction", .hasArg = true,.foundArg = &phiArg,	.foundOpt = NULL },
		{ .opt = "R", .optAlternative = "mean-particle-radius", .hasArg = true,.foundArg = &RsArg,	.foundOpt = NULL },
		{ .opt = "G", .optAlternative = "shear-modulus", .hasArg = true,.foundArg = &GArg,		.foundOpt = NULL },
		{ .opt = "B", .optAlternative = "burgers-vector", .hasArg = true,.foundArg = &bArg,		.foundOpt = NULL },
		{ .opt = "m", .optAlternative = "taylor-factor", .hasArg = true,.foundArg = &MArg,		.foundOpt = NULL },
		{0},
	};

	if (parseArgs(argc, argv, &arguments->common, options) != kCommonConstantReturnTypeSuccess)
	{
		fprintf(stderr, "Parsing command line arguments failed\n");
		printUsage();

		return kCommonConstantReturnTypeError;
	}

	if (arguments->common.isHelpEnabled)
	{
		printUsage();

		exit(EXIT_SUCCESS);
	}

	if (arguments->common.isOutputSelected)
	{
		fprintf(stderr, "Error: Output select option not supported.\n");

		return kCommonConstantReturnTypeError;
	}

	if (arguments->common.isInputFromFileEnabled && arguments->common.isMonteCarloMode)
	{
		fprintf(stderr, "Error: Reading from an input file is not supported for Monte Carlo mode.\n");

		return kCommonConstantReturnTypeError;
	}

	if (gammaArg != NULL)
	{
		double gamma;

		if (arguments->common.isMonteCarloMode)
		{
			if (strstr(gammaArg, kConstantStringUx) != NULL)
			{
				fprintf(stderr, "Error: Native Monte Carlo is not compatible with Ux strings from command line.\n");

				return kCommonConstantReturnTypeError;
			}
		}

		int ret = parseDoubleChecked(gammaArg, &gamma);

		if (ret != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The gamma must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->gamma = gamma;
	}

	if (phiArg != NULL)
	{
		double phi;

		if (arguments->common.isMonteCarloMode)
		{
			if (strstr(phiArg, kConstantStringUx) != NULL)
			{
				fprintf(stderr, "Error: Native Monte Carlo is not compatible with Ux strings from command line.\n");

				return kCommonConstantReturnTypeError;
			}
		}

		int ret = parseDoubleChecked(phiArg, &phi);

		if (ret != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The phi must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->phi = phi;
	}

	if (RsArg != NULL)
	{
		double Rs;

		if (arguments->common.isMonteCarloMode)
		{
			if (strstr(RsArg, kConstantStringUx) != NULL)
			{
				fprintf(stderr, "Error: Native Monte Carlo is not compatible with Ux strings from command line.\n");

				return kCommonConstantReturnTypeError;
			}
		}

		int ret = parseDoubleChecked(RsArg, &Rs);

		if (ret != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The Rs must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->Rs = Rs;
	}

	if (GArg != NULL)
	{
		double G;

		if (arguments->common.isMonteCarloMode)
		{
			if (strstr(GArg, kConstantStringUx) != NULL)
			{
				fprintf(stderr, "Error: Native Monte Carlo is not compatible with Ux strings from command line.\n");

				return kCommonConstantReturnTypeError;
			}
		}

		int ret = parseDoubleChecked(GArg, &G);

		if (ret != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The G must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->G = G;
	}

	if (bArg != NULL)
	{
		double b;

		if (arguments->common.isMonteCarloMode)
		{
			if (strstr(bArg, kConstantStringUx) != NULL)
			{
				fprintf(stderr, "Error: Native Monte Carlo is not compatible with Ux strings from command line.\n");

				return kCommonConstantReturnTypeError;
			}
		}

		int ret = parseDoubleChecked(bArg, &b);

		if (ret != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The b must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->b = b;
	}

	if (MArg != NULL)
	{
		double M;

		if (arguments->common.isMonteCarloMode)
		{
			if (strstr(MArg, kConstantStringUx) != NULL)
			{
				fprintf(stderr, "Error: Native Monte Carlo is not compatible with Ux strings from command line.\n");

				return kCommonConstantReturnTypeError;
			}
		}

		int ret = parseDoubleChecked(MArg, &M);

		if (ret != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The M must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->M = M;
	}

	return kCommonConstantReturnTypeSuccess;
}

void
loadInputs(
	double *		gamma,
	double *		phi,
	double *		Rs,
	double *		G,
	double *		b,
	double *		M,
	double *		inputDistributions,
	CommandLineArguments *	arguments)
{
	if (arguments->common.isInputFromFileEnabled)
	{
		*b	= inputDistributions[kInputDistributionIndexB];
		*G	= inputDistributions[kInputDistributionIndexG];
		*gamma	= inputDistributions[kInputDistributionIndexGamma];
		*M	= inputDistributions[kInputDistributionIndexM];
		*phi	= inputDistributions[kInputDistributionIndexPhi];
		*Rs	= inputDistributions[kInputDistributionIndexRs];
	}
	else
	{
		*b	= arguments->b;
		*G	= arguments->G;
		*gamma	= arguments->gamma;
		*M	= arguments->M;
		*phi	= arguments->phi;
		*Rs	= arguments->Rs;
	}

	return;
}

void
printJSONFormattedOutput(
	double			sigmaCMpa,
	double			cpuTimeUsedInSeconds,
	CommandLineArguments *	arguments)
{
	if (arguments->common.isTimingEnabled)
	{
		JSONVariable variables[2] = {
			{
				.variableSymbol = "sigmaCMpa",
				.variableDescription = "Cutting stress (σc)",
				.values = (JSONVariablePointer) { .asDouble = &sigmaCMpa},
				.type = kJSONVariableTypeDouble,
				.size = 1,
			},
			{
				.variableSymbol = "cpuTimeUsed",
				.variableDescription = "CPU time used (s)",
				.values = (JSONVariablePointer) { .asDouble = &cpuTimeUsedInSeconds},
				.type = kJSONVariableTypeDoubleParticle,
				.size = 1,
			}
		};

		printJSONVariables(variables, 2, "Precipitate \\\"cutting\\\" dislocation model from Brown and Ham");
	}
	else
	{
		JSONVariable variables[1] = {
			{
				.variableSymbol = "sigmaCMpa",
				.variableDescription = "Cutting stress (σc)",
				.values = (JSONVariablePointer) { .asDouble = &sigmaCMpa},
				.type = kJSONVariableTypeDouble,
				.size = 1,
			}
		};
		printJSONVariables(variables, 1, "Precipitate \\\"cutting\\\" dislocation model from Brown and Ham");
	}

	return;
}
