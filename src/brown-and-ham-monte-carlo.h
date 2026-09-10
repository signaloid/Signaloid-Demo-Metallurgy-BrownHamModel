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
 *	@brief	Calculate and return the cutting stress for each Monte Carlo
 *		sample, with one independent draw of the model's inputs per
 *		iteration, and write the result into `monteCarloOutputSamples`.
 *		Each sample is produced by `brownHamModelMonteCarloSample()`
 *		(declared in `kernel.h`).
 *
 *	@param	arguments			: Command-line arguments, including the override flags.
 *	@param	numberOfMonteCarloIterations	: Number of Monte Carlo iterations to run.
 *	@param	monteCarloOutputSamples		: Array of `numberOfMonteCarloIterations` doubles to fill with samples.
 *	@return	double				: Returns the last element of `monteCarloOutputSamples`.
 */
double
computeSigmaMonteCarlo(
	CommandLineArguments *  arguments,
	size_t                  numberOfMonteCarloIterations,
	double *                monteCarloOutputSamples);
