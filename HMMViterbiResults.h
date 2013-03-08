/*
 * HMMViterbiResults.h
 *
 *	This is the header file for the HMMViterbiResults object. 
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

#ifndef HMMVITERBIRESULTS_H
#define HMMVITERBIRESULTS_H
#include "HMMProbabilities.h"
#include <map>
#include <vector>
#include <string>
using namespace std;

class HMMViterbiResults
{
public:
	// Constuctors
	// ==============================================
	HMMViterbiResults();
	HMMViterbiResults(int iteration, int numberOfStates);

	// Destructor
	// =============================================
	~HMMViterbiResults();

	// Public Attributes
	int iteration;
	int numStates;
	vector<int> stateCounts;
	vector<int> segmentCounts;
	map<int,vector<pair<int,int>>> segments;
	map<int, map<string, int>> emissionCounts;
	vector<vector<int>> transitionCounts;
	HMMProbabilities* probabilities;

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
	string resultsWithoutSegments();

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
	string allResults();

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
	void calculateProbabilities(HMMProbabilities* previousProbs);

private:

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
	string stateHistogramResultsString();

	// string segmentHistogramResultsString()
	//  Purpose:
	//		Returns a string representing the segment histogram
	//
	//		format:
	//			<result type="segment_histogram">
	//				<<state>>=<<segment count>>,
	//			</result>
	string segmentHistogramResultsString();

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
	string probabilitiesResultsString();

	// string segmentResultsString()
	//  Purpose:
	//		Returns a string representing the segments
	//
	//		format:
	//			<result type="segments">
	//				(segment1start, segment1end),(segment2start, segment2end),...
	//			</result>
	string segmentResultsString();

	// string statesResultsString()
	//  Purpose:
	//		Returns a string representing the states
	//
	//		format:
	//			<result type="states">
	//				<<state1>>,<<state2>>,...
	//			</result>
	string statesResultsString();

	// string intitiationProbabiltiesResultsString()
	//  Purpose:
	//		Returns a string representing the initiation probablities
	//
	//		format:
	//			<result type="initiation_probabilites">
	//				<<state>>=<<initiation probability>>,
	//			</result>
	string intitiationProbabiltiesResultsString();

	// transitionProbablitiesResultsString(int state)
	//  Purpose:
	//		Returns a string representing the transition probablities for a state
	//
	//		format:
	//			<result type="transition_probabilites" state="<<state>>">
	//				<<to state>>=<<transition probability>>,
	//			</result>
	string transitionProbablitiesResultsString(int state);

	// emissionProbablitiesResultsString(int state)
	//  Purpose:
	//		Returns a string representing the emission probablities for a state
	//
	//		format:
	//			<result type="emission_probabilites" state="<<state>>">
	//				<<residue>>=<<emission probability>>,
	//			</result>
	string emissionProbablitiesResultsString(int state);

	// transitionCountsResultsString(int state)
	//  Purpose:
	//		Returns a string representing the transimission counts
	//
	//		format:
	//			<result type="transitionCounts">
	//				<<transition>>=<<transition count>>,
	//			</result>
	string transitionCountsResultsString();
	
	struct sort_segment {
		bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
			return left.second - left.first < right.second - right.first;
		}
	};

};

#endif // HMMVITERBIRESULTS_H
