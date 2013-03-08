/*
 * HMMViterbiResults.cpp
 *
 *	This is the cpp file for the HMMViterbiResults object. 
 *  HMMViterbiResults contains the rsults for one iteration of viterbi
 *  training for a hidden markov model.
 *
 *  Important Attributes
 *
 *		stateCounts - how many times a state occurs in the viterbi path
 *		segmentCounts - how many segments (i.e., continuos occurencee of
 *						one state) of each state occurs
 *		segments - collection of start and stop values for each segment
 *		transitionCounts - counts for how many time states transition (both
 *						   from one state to another and from one state to 
 *						   the stame state)
 *		probabilities - the probabilites that are calculated from the above 
 *						results
 *
 *	resultsWithoutSegments() - returns a string of all results except segments
 *	allResults() - returns a string of all results (including segements)
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */
#include "HMMViterbiResults.h"
#include "StringUtilities.h"
#include <sstream>
#include <vector>
#include <algorithm>

// Constuctors
// ==============================================
HMMViterbiResults::HMMViterbiResults() {
}

HMMViterbiResults::HMMViterbiResults(int anIteration, int numberOfStates) {

	iteration = anIteration;
	numStates = numberOfStates;
	probabilities = new HMMProbabilities(numStates);


	// initialize counts vectors
	for (int i = 0; i < numStates; i++) {
		stateCounts.push_back(0);
		segmentCounts.push_back(0);
		segments[i] = vector<pair<int,int>>();
		transitionCounts.push_back(vector<int>());
		for (int j = 0; j < numStates; j++) {
			transitionCounts[i].push_back(0);
		}
		if (i > 0) {
			for (pair<string, int> mapPair : probabilities->emissionResidueMap) {
				string& residue = mapPair.first;
				emissionCounts[i][residue] = 0;
			}
		}

	}
}

// Destructor
// =============================================
HMMViterbiResults::~HMMViterbiResults(){
}

// Public Methods
// =============================================

// string resultsWithoutSegments()
//  Purpose:
//		Returns a string representing the viterbi results for 
//		a particular iteration
//
//		format:
//			<result type="viterbi_iteration" iteration="<< iteration >>">
//				<<stateHistogramResultsString>>
//				<<segmentHistogramResultsString>>
//				<<probabiltiesResultsString>>
//			</result>
string HMMViterbiResults::resultsWithoutSegments() {
	stringstream ss;

	// Header
	ss << "    <result type=\"viterbi_iteration\" iteration=\"" << iteration << "\">\n";

	// Results
	ss 
		<< stateHistogramResultsString()
		<< segmentHistogramResultsString()
		<< probabilitiesResultsString();

//	ss 	<< transitionCountsResultsString();

	// Footer
	ss << "    </result>\n";

	return ss.str();
}

// string allResults()
//  Purpose:
//		Returns a string representing the viterbi results for 
//		a particular iteration
//
//		format:
//			<result type="viterbi_iteration" iteration="<< iteration >>">
//				<<stateHistogramResultsString>>
//				<<segmentHistogramResultsString>>
//				<<probabiltiesResultsString>>
//			</result>
//			<<segmentResultsString>>
string HMMViterbiResults::allResults() {
	stringstream ss;

	// Results
	ss 
		<< resultsWithoutSegments()
		<< segmentResultsString();

	return ss.str();
}

// calculateProbabilities(HMMProbabilities* previousProbs)
//  Purpose:
//		Calculates the probabilties using the information stored in
//		this object as well as the probabilites from the previous
//		iteration of the viterbi training.
//
//		For now, only the transistion probabilties are being recalculated
//		from the data in this object.  The initaition and emission
//		probabilties are being held steady and are set to the values
//		from the previous iteration
//	Postconditions:
//		probabilites - will be populated
void HMMViterbiResults::calculateProbabilities(HMMProbabilities* previousProbs) {

	// initiation probabilties - Use initation from previous probabilites
	for (int state = 1; state < numStates; state++) {
		probabilities->setInitiationProbability(state, previousProbs->initiationProbability(state));
	}

	// emission probabilities
	for (int state = 1; state < numStates; state++) {
		for (pair<string, int> mapPair : probabilities->emissionResidueMap) {
			string& residue = mapPair.first;
			long double newProbability = 
				emissionCounts.at(state).at(residue) / (double) stateCounts[state];
			probabilities->setEmissionProbability(state, residue, newProbability);
		}
	}

	// transition probabilites
	for (int firstState = 1; firstState < numStates; firstState++) {
		for (int secondState = 1; secondState < numStates; secondState++) {
			long double newProbability =
				transitionCounts[firstState][secondState] / (double) stateCounts[firstState];
			probabilities->setTransitionProbability(firstState, secondState, newProbability);
		}
	}
}

// Private Methods
// =============================================

// string stateHistogramResultsString()
//  Purpose:
//		Returns a string representing the state histogram
//
//		format:
//			<result type="state_histogram">
//				<<state>>=<<state count>>,
//			</result>
string HMMViterbiResults::stateHistogramResultsString() {
	stringstream ss;

	for (int i = 0; i < numStates; i++) {
		ss 
			<< i + 1 
			<< "="
			<< stateCounts[i];

		if ( i < numStates -1)
		   ss << ",";
	}

	return StringUtilities::xmlResult("state_histogram", ss.str());
}

// string segmentHistogramResultsString()
//  Purpose:
//		Returns a string representing the segment histogram
//
//		format:
//			<result type="segment_histogram">
//				<<state>>=<<segment count>>,
//			</result>
string HMMViterbiResults::segmentHistogramResultsString() {
	stringstream ss;

	for (int i = 0; i < numStates; i++) {
		ss 
			<< i + 1 
			<< "="
			<< segmentCounts[i];

		if ( i < numStates -1)
		   ss << ",";
	}

	return StringUtilities::xmlResult("segment_histogram", ss.str());
}

// string probabilitiesResultsString()
//  Purpose:
//		Returns a string representing the probabilites
//
//		format:
//			<<statesResultsString>>
//			<<initiationProbabilitesResultsString>>
//			<<transmissionProbabilitesResultsString>>
//			...
//			<<emissionProbabilitesResultsString>>
//			...
string HMMViterbiResults::probabilitiesResultsString() {
	return probabilities->probabilitiesResultsString();
}

// string segmentResultsString()
//  Purpose:
//		Returns a string representing the segments
//
//		format:
//			<result type="segments">
//				(segment1start, segment1end),(segment2start, segment2end),...
//			</result>
string HMMViterbiResults::segmentResultsString() {
	stringstream ss;

	vector<pair<int, int>>& gcSegments = segments[2];
	sort(gcSegments.begin(), gcSegments.end(), sort_segment());
	
	int counter = 0;
	for (int i = gcSegments.size() - 1; i >= 0; i--) {
		pair<int, int>& segment = gcSegments[i];
		ss 
			<< "("
			<< segment.first
			<< ","
			<< segment.second
			<< ")";
		if (i > 0)
			ss << ",";

		counter++;
		if (counter % 5 == 0)
			ss << "\n";
	}

	return StringUtilities::xmlResult("segment_list", ss.str());
}

// transitionCountsResultsString(int state)
//  Purpose:
//		Returns a string representing the transimission counts
//
//		format:
//			<result type="transitionCounts">
//				<<transition>>=<<transition count>>,
//			</result>
string HMMViterbiResults::transitionCountsResultsString() {
	stringstream ss;

	// Header 
	ss << "        <transition_counts>";

	// States
	for (int i = 0; i < numStates; i++) {
		for (int j = 0; j < numStates; j++) {
			ss
				<< i + 1 << j + 1
				<< "="
				<< transitionCounts[i][j];

			if ( i < numStates - 1 || j < numStates - 1)
			   ss << ",";
		}
	}

	// Footer
	ss << "</transition_counts>\n";

	return ss.str();
}
