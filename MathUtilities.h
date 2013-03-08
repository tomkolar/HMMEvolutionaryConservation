/*
 * MathUtilieis.h
 *
 *  The MathUtilities object is a container for extended mathmatcial operations.
 *  Initially it contains extensions for logarithm functions so that a zero
 *  input can be handled.
 *
 *	The logarithm extensions are based on the paper "Numerically Stable Hidden
 *  Markov Model Implementation" by Tobias Mann. 
 *		//bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
 *
 *
 *  Created on: 2-22-13
 *      Author: tomkolar
 */

#ifndef MATHUTILITITES_H
#define STRINGUTMATHUTILITITES_HILITIES_H_

using namespace std;

class MathUtilities
{
public:

	// Constuctors
	// ==============================================
	MathUtilities();

	// Destructor
	// =============================================
	~MathUtilities();

	// Public Class Methods
	// =============================================

	// double eexp (double aNumber)
	//	Purpose:
	//		Extends the standard exponential function exp() to allow zero as
	//		an input.
	static long double eexp (long double aNumber);

	// double eln (double aNumber)
	//	Purpose:
	//		Extends the standar logarithm function ln() to allow zero as
	//		an input.
	static long double eln (long double aNumber);

	// double elnsum (double lnOfX, double lnOfY)
	//	Purpose:
	//		Returns the ln of the sum of two ln values allowing lnOfX or
	//		lnOfY to be zero
	static long double elnsum (long double lnOfX, long double lnOfY);

	// double elnprod (double lnOfX, double lnOfY)
	//	Purpose:
	//		Returns the ln of the product of two ln values allowing lnOfX or
	//		lnOfY to be zero
	static long double elnprod (long double lnOfX, long double lnOfY);

	// bool isNaN(double aValue)
	//	Purpose:  Returns true if aValue is not a number
	static bool isNaN(long double var);
};

#endif // MATHUTILITITES_H
