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

#pragma once

#include <stdlib.h>
#include "utilities.h"

/**
 *	@brief	Evaluate the Brown and Ham precipitate "cutting" dislocation model
 *		formula for one realization of its six inputs and return the
 *		predicted cutting stress. The formula is pure arithmetic (`sqrt`,
 *		`pow`, and the basic operators), so it is shared unmodified between
 *		UxHw mode (where `gamma`, `phi`, `Rs`, `G`, and `M` may each carry a
 *		full probability distribution) and Monte Carlo mode (where each is
 *		a single sample). No UxHw API calls are required here: Signaloid
 *		hardware propagates distributions through plain arithmetic without
 *		explicit distributional API calls.
 *
 *		                    ⎛    _________________    ⎞
 *		       ⎛ M ⋅ γ  ⎞   ⎜   ╱8.0 ⋅ γ ⋅ φ ⋅ Rs     ⎟
 *		  σ  = ⎜─────── ⎟ ⋅ ⎜  ╱ ───────────────── - φ⎟
 *		   c   ⎝2.0 ⋅ b ⎠   ⎝╲╱  π ⋅ G ⋅ pow(b, 2)    ⎠
 *
 *	@param	gamma	: Anti-phase boundary energy (APB energy), units J/m^2.
 *	@param	phi	: Precipitate volume fraction.
 *	@param	Rs	: Mean particle radius on plane, units m.
 *	@param	G	: Shear modulus, units Pa.
 *	@param	b	: Magnitude of the Burgers vector, units m.
 *	@param	M	: Taylor factor.
 *	@return	double	: Returns the predicted cutting stress, `σc`.
 */
double
computeBrownHamModelOutput(
	double  gamma,
	double  phi,
	double  Rs,
	double  G,
	double  b,
	double  M);

/**
 *	@brief	Draw one independent Monte Carlo sample of the model's inputs and
 *		return the resulting cutting stress.
 *
 *	@param	arguments	: Command-line arguments, including the override flags.
 *	@return	double		: Returns the cutting stress for this sample.
 */
double
brownHamModelMonteCarloSample(CommandLineArguments * arguments);

/**
 *	@brief	UxHw-mode calculation kernel. Computes the cutting stress using
 *		distributional arithmetic on the single already-resolved input set
 *		in `arguments` (populated from the command line, its defaults, or
 *		an input CSV file by `getCommandLineArguments()` / `loadInputs()`).
 *		Writes the result into `outputVariables[kOutputDistributionIndexSigma]`
 *		and into `monteCarloOutputSamples[0]`.
 *
 *	@param	arguments		: Command-line arguments.
 *	@param	outputVariables		: Array of size `kOutputDistributionIndexMax` to fill.
 *	@param	monteCarloOutputSamples	: Single-element array for the distributional result.
 *	@return	double			: Returns the cutting stress, `σc`.
 */
double
calculateOutputUxHw(
	CommandLineArguments *  arguments,
	double *                outputVariables,
	double *                monteCarloOutputSamples);

/**
 *	@brief	Monte Carlo calculation kernel. Runs
 *		`arguments->common.numberOfMonteCarloIterations` independent samples
 *		into `monteCarloOutputSamples`, using `brownHamModelMonteCarloSample`
 *		for each one (no UxHw distributional API calls in this file itself).
 *		Writes the last sample's value into
 *		`outputVariables[kOutputDistributionIndexSigma]`.
 *
 *	@param	arguments		: Command-line arguments.
 *	@param	outputVariables		: Array of size `kOutputDistributionIndexMax` to fill.
 *	@param	monteCarloOutputSamples	: Array of `numberOfMonteCarloIterations` doubles, filled with samples.
 *	@return	double			: Returns the cutting stress for the last sample.
 */
double
calculateOutputMonteCarlo(
	CommandLineArguments *  arguments,
	double *                outputVariables,
	double *                monteCarloOutputSamples);
