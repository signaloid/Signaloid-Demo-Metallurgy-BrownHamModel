/*
 *	Copyright (c) 2020--2022, Signaloid.
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
#include <stdint.h>
#include <stdio.h>
#include <uxhw.h>

/*
 *	Precipitate "cutting" dislocation model from Brown and Ham
 *
 *	Overview:
 *	Models an equation for a materials precipitate "cutting" dislocation model from Brown and Ham.
 *
 *	Inputs:
 *	The inputs and their ranges are:
 *	-	`gamma`:	0.15 to 0.25
 *	-	`phi`:		0.30 to 0.45
 *	-	`Rs`:		1 10^-8 to 3 10^-8
 *	-	`G`:		6 10^10 to 8 10^10
 *	-	`b`:		2.54 10^-9 to 2.54 10^-9 (i.e., constant)
 *	-	`M`:		2.9 to 3.2.
 *
 *	The parameter `gamma` is the APB energy with units J/m^2, `phi` is the 
 *	precipitate volume fraction, `Rs` is mean particle radius on plane with units m,
 *	`G` is the shear modulus with units Pa, `b` is the magnitude of the Burgers
 *	vector with units m, and `M` is the Taylor factor.
 *
 *	Outputs:
 *	The output is the cutting stress (alloy strength), `σc` where
 *
 *                    ⎛    _________________    ⎞
 *       ⎛ M ⋅ γ  ⎞   ⎜   ╱8.0 ⋅ γ ⋅ φ ⋅ Rs     ⎟
 *  σ  = ⎜─────── ⎟ ⋅ ⎜  ╱ ───────────────── - φ⎟
 *   c   ⎝2.0 ⋅ b ⎠   ⎝╲╱  π ⋅ G ⋅ pow(b, 2)    ⎠
 *
 */

static void
loadInputs(double *  G, double *  M, double *  Rs, double *  b, double *  gamma, double *  phi)
{
	double empiricalTaylorFactorValues[] = {
						3.2,
						3.9,
						4.1,
						3.2,
						3.8,
						3.8,
						2.1,
						3.0,
						1.9,
						3.9,
						2.3,
						2.2,
						3.2,
						2.2,
						3.9,
						2.2,
						1.9,
						3.2,
						3.9,
						3.1,
					};

	for (int i = 0; i < sizeof(empiricalTaylorFactorValues)/sizeof(double); i++)
	{
		*M += empiricalTaylorFactorValues[i];
	}
	*M /= sizeof(empiricalTaylorFactorValues)/sizeof(double);

	*G		= 7E10;
	*Rs		= 2E-8;
	*b		= 2.54E-10;
	*gamma	= 0.2;
	*phi	= 0.375;
}

int
main(int argc, char *	argv[])
{
	double	G, M, Rs, b, gamma, phi, sigmaCMpa;


	loadInputs(&G, &M, &Rs, &b, &gamma, &phi);

	/*
	 *                    ⎛    _________________    ⎞
	 *       ⎛ M ⋅ γ  ⎞   ⎜   ╱8.0 ⋅ γ ⋅ φ ⋅ Rs     ⎟
	 *  σ  = ⎜─────── ⎟ ⋅ ⎜  ╱ ───────────────── - φ⎟
	 *   c   ⎝2.0 ⋅ b ⎠   ⎝╲╱  π ⋅ G ⋅ pow(b, 2)    ⎠
	 */
	sigmaCMpa = ((M*gamma)/(2.0*b))*(sqrt((8.0*gamma*phi*Rs)/(M_PI*G*pow(b, 2))) - phi)/1000000;

	printf("Alloy strength (σc)\t\t= %.1E MPa\n", sigmaCMpa);


	return 0;
}
