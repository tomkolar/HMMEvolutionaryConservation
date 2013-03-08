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
#include "MathUtilities.h"
#include <math.h>
#include <stdexcept>
#include <limits>

// Constuctors
// ==============================================
MathUtilities::MathUtilities() {
}

// Destructor
// =============================================
MathUtilities::~MathUtilities() {
}

// Public Class Methods
// =============================================

// double eexp (double aNumber)
//	Purpose:
//		Extends the standard exponential function exp() to allow zero as
//		an input.
long double MathUtilities::eexp (long double aNumber) {

	if (isNaN(aNumber)) 
		return 0;

	return exp(aNumber);
}

// double eln (double aNumber)
//	Purpose:
//		Extends the standar logarithm function ln() to allow zero as
//		an input.
long double MathUtilities::eln (long double aNumber) {

	if (aNumber == 0.0) 
		return std::numeric_limits<double>::quiet_NaN();

	if (aNumber > 0.0)
		return log(aNumber);

	throw out_of_range("Negative Input Error");
}

// double elnsum (double lnOfX, double lnOfY)
//	Purpose:
//		Returns the ln of the sum of two ln values allowing lnOfX or
//		lnOfY to be zero
long double MathUtilities::elnsum (long double lnOfX, long double lnOfY) {
	if (isNaN(lnOfX) || isNaN(lnOfY)) {
		if (isNaN(lnOfX))
			return lnOfY;
		else
			return lnOfX;
	}

	if (lnOfX > lnOfY) 
		return lnOfX + eln(1.0 + exp(lnOfY - lnOfX));

	return lnOfY + eln(1.0 + exp(lnOfX - lnOfY));
}

// double elnprod (double lnOfX, double lnOfY)
//	Purpose:
//		Returns the ln of the product of two ln values allowing lnOfX or
//		lnOfY to be zero
long double MathUtilities::elnprod (long double lnOfX, long double lnOfY) {

	if (isNaN(lnOfX))
		return lnOfX; // return NAN

	if (isNaN(lnOfY))
		return lnOfY; // return NAN

	return lnOfX + lnOfY;
}

// bool isNaN(double aValue)
//	Purpose:  Returns true if aValue is not a number
bool MathUtilities::isNaN(long double aValue)
{
	volatile long double d = aValue;
	return d != d;
}
