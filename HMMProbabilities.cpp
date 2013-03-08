/*
 * HMMProbablities.cpp
 *
 *	This is the cpp file for the HMMProbabilties object. HMMProbabilities
 *  is a collection of all probabilties needed for a hidden markov model. This
 *  includes initiation, emission and transition probabilties.  There are 
 *	convenience methods for setting and retriving probabilties as well as
 *  the log value of each probabilty.
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */
#include "HMMProbabilities.h"
#include "StringUtilities.h"
#include <cmath>
#include <sstream>
#include <string>
#include <limits>
#include <iostream>
#include <fstream>

// Constuctors
// ==============================================
HMMProbabilities::HMMProbabilities() {
}

HMMProbabilities::HMMProbabilities(int numOfStates) {
	numStates = numOfStates;
	createEmissionResidueMap();

	// Initialize all probabilities to zero
	for (int i = 0; i < numStates; i++) {
		setInitiationProbability(i, 0);
		for (int j = 0; j < numStates; j++) {
			setTransitionProbability(i, j, 0);
		}
		for (pair<string, int> mapPair : emissionResidueMap) {
			string& residue = mapPair.first;
			setEmissionProbability(i, residue, 0);
		}
	}
}


// Destructor
// =============================================
HMMProbabilities::~HMMProbabilities() {
}

// Public Class Methods
// =============================================

// HMMProbabilities* initialProbabilities()
//  Purpose: 
//		Returns a probabilites object initialzed to the initial
//		probabilites required by genome540 homework #7	
HMMProbabilities* HMMProbabilities::initialProbabilities(
	string neutralCountsFile, string conservedCountsFile) {

	HMMProbabilities* probs = new HMMProbabilities(3);

	// initiation probabilties
	probs->setInitiationProbability(1, 0.95);
	probs->setInitiationProbability(2, 0.05);

	// transition probabilities
	probs->setTransitionProbability(1, 1, 0.95);
	probs->setTransitionProbability(1, 2, 0.05);
	probs->setTransitionProbability(2, 1, 0.10);
	probs->setTransitionProbability(2, 2, 0.90);

	// emission probabilities
	probs->populateEmissionProbabilities(1, neutralCountsFile);
	probs->populateEmissionProbabilities(2, conservedCountsFile);

	return probs;

}

// Public Methods
// =============================================

// double emissionProbability(int state, char residue)
//  Purpose: 
//		Returns the emission probability for the state and residue
long double HMMProbabilities::emissionProbability(int state, string residue) {
	return emissionProbabilities.at(state).at(getEmissionResidueIndex(residue));
}

// double initiationProbability(int state)
//  Purpose: 
//		Returns the initiation probability for the state
long double HMMProbabilities::initiationProbability(int state) {
	return initiationProbabilities[state];
}
	
// double transitionProbability(int beginState, int endState)
//  Purpose: 
//		Returns the transition probability for transition from beginState
//		to endState
long double HMMProbabilities::transitionProbability(int beginState, int endState) {
	return transitionProbabilities[beginState][endState];
}

// double logEmissionProbability(int state, char residue)
//  Purpose: 
//		Returns the log of the emission probability for the state and residue
long double HMMProbabilities::logEmissionProbability(int state, string residue) {
	return logEmissionProbabilities.at(state).at(getEmissionResidueIndex(residue));
}

// double logInitiationProbability(int state)
//  Purpose: 
//		Returns the log of the initiation probability for the state
long double HMMProbabilities::logInitiationProbability(int state) {
	return logInitiationProbabilities[state];
}
	
// double logTransitionProbability(int beginState, int endState)
//  Purpose: 
//		Returns the log of the transition probability for transition from
//		beginState to endState
long double HMMProbabilities::logTransitionProbability(int beginState, int endState) {
	return logTransitionProbabilities[beginState][endState];
}

// setEmissionProbability(int state, char residue, double value)
//  Purpose: 
//		Sets the emission probability for the state and residue to value
//	Postconditions:
//		emissionProbabilites - value set for state/residue
//		logEmissionProbabilites - value set for state/residue
void HMMProbabilities::setEmissionProbability(int state, string residue, long double value) {
	emissionProbabilities[state][getEmissionResidueIndex(residue)] = value;
	double logVal;
	if (value == 0)
		logVal = std::numeric_limits<double>::quiet_NaN();
	else
		logVal = log(value);
	logEmissionProbabilities[state][getEmissionResidueIndex(residue)] = logVal;
}

// setInitiationProbability(int state, double value)
//  Purpose: 
//		Sets the initiation probability for the state to value
//	Postconditions:
//		initiationProbabilites - value set for state
//		logInitiationProbabilites - value set for state
void HMMProbabilities::setInitiationProbability(int state, long double value) {
	initiationProbabilities[state] = value;
	double logVal;
	if (value == 0)
		logVal = std::numeric_limits<double>::quiet_NaN();
	else
		logVal = log(value);
	logInitiationProbabilities[state] = logVal;
}

// setTransitionProbability(int beginState, int endState, double value)
//  Purpose: 
//		Sets the transition probability from the beginState to the endState
//		to value
//	Postconditions:
//		transitionProbabilites - value set for beginState to endState
//		logTransitionProbabilites - value set for beginState to endState
void HMMProbabilities::setTransitionProbability(int beginState, int endState, long double value) {
	transitionProbabilities[beginState][endState] = value;
	double logVal;
	if (value == 0)
		logVal = std::numeric_limits<double>::quiet_NaN();
	else
		logVal = log(value);
	logTransitionProbabilities[beginState][endState] = logVal;
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
string HMMProbabilities::probabilitiesResultsString() {
	stringstream ss;

	// Begin Model
	ss << "      <model type=\"hmm\">\n";

	// States
	ss << statesResultsString();

	// Probabiltiies
	ss << intitiationProbabiltiesResultsString();
	for (int i = 1; i < numStates; i++)
		ss << transitionProbablitiesResultsString(i);
	for (int i = 1; i < numStates; i++)
		ss << emissionProbablitiesResultsString(i);
	
	// End Model
	ss << "      </model>\n";

	return ss.str();
}

// string statesResultsString()
//  Purpose:
//		Returns a string representing the states
//
//		format:
//			<result type="states">
//				<<state1>>,<<state2>>,...
//			</result>
string HMMProbabilities::statesResultsString() {
	stringstream ss;

	// Header 
	ss << "        <states>";

	// States
	for (int i = 1; i < numStates; i++) {
		ss << i ;

		if ( i < numStates - 1)
		   ss << ",";
	}

	// Footer
	ss << "</states>\n";

	return ss.str();
}

// string intitiationProbabiltiesResultsString()
//  Purpose:
//		Returns a string representing the initiation probablities
//
//		format:
//			<result type="initiation_probabilites">
//				<<state>>=<<initiation probability>>,
//			</result>
string HMMProbabilities::intitiationProbabiltiesResultsString() {
	stringstream ss;
	ss.precision(5);

	// Header 
	ss << "        <initial_state_probabilities>";

	// States
	for (int i = 1; i < numStates; i++) {
		ss
			<< i 
			<< "="
			<< initiationProbability(i);

		if ( i < numStates - 1)
		   ss << ",";
	}

	// Footer
	ss << "</initial_state_probabilities>\n";

	return ss.str();
}

// transitionProbablitiesResultsString(int state)
//  Purpose:
//		Returns a string representing the transition probablities for a state
//
//		format:
//			<result type="transition_probabilites" state="<<state>>">
//				<<to state>>=<<transition probability>>,
//			</result>
string HMMProbabilities::transitionProbablitiesResultsString(int state) {
	stringstream ss;
	ss.precision(5);

	// Header 
	ss << "        <transition_probabilities state=\"" << state << "\">";

	// States
	for (int i = 1; i < numStates; i++) {
		ss
			<< i
			<< "="
			<< transitionProbability(state, i);

		if ( i < numStates - 1)
		   ss << ",";
	}

	// Footer
	ss << "</transition_probabilities>\n";

	return ss.str();
}

// emissionProbablitiesResultsString(int state)
//  Purpose:
//		Returns a string representing the emission probablities for a state
//
//		format:
//			<result type="emission_probabilites" state="<<state>>">
//				<<residue>>=<<emission probability>>,
//			</result>
string HMMProbabilities::emissionProbablitiesResultsString(int state) {
	stringstream ss;
	ss.precision(5);

	// Header 
	ss << "        <emission_probabilities state=\"" << state << "\">";

	// Residues
	for (pair<string, int> mapPair : emissionResidueMap) {
		string& residue = mapPair.first;
		ss << residue << "=" << emissionProbability(state, residue) << ",";
	}

	// Footer
	ss << "</emission_probabilities>\n";

	return ss.str();
}

// map<string, int> createEmissionMap()
//  Purpose: 
//		Creates a map of the index location for a trinucleotide emission
//		in the emission probabilities array
void HMMProbabilities::createEmissionResidueMap() {

	emissionResidueMap["AAA"]= 0;
	emissionResidueMap["AAC"]= 1;
	emissionResidueMap["AAT"]= 2;
	emissionResidueMap["AAG"]= 3;
	emissionResidueMap["AA-"]= 4;
	emissionResidueMap["ACA"]= 5;
	emissionResidueMap["ACC"]= 6;
	emissionResidueMap["ACT"]= 7;
	emissionResidueMap["ACG"]= 8;
	emissionResidueMap["AC-"]= 9;
	emissionResidueMap["ATA"]= 10;
	emissionResidueMap["ATC"]= 11;
	emissionResidueMap["ATT"]= 12;
	emissionResidueMap["ATG"]= 13;
	emissionResidueMap["AT-"]= 14;
	emissionResidueMap["AGA"]= 15;
	emissionResidueMap["AGC"]= 16;
	emissionResidueMap["AGT"]= 17;
	emissionResidueMap["AGG"]= 18;
	emissionResidueMap["AG-"]= 19;
	emissionResidueMap["A-A"]= 20;
	emissionResidueMap["A-C"]= 21;
	emissionResidueMap["A-T"]= 22;
	emissionResidueMap["A-G"]= 23;
	emissionResidueMap["A--"]= 24;
	emissionResidueMap["CAA"]= 25;
	emissionResidueMap["CAC"]= 26;
	emissionResidueMap["CAT"]= 27;
	emissionResidueMap["CAG"]= 28;
	emissionResidueMap["CA-"]= 29;
	emissionResidueMap["CCA"]= 30;
	emissionResidueMap["CCC"]= 31;
	emissionResidueMap["CCT"]= 32;
	emissionResidueMap["CCG"]= 33;
	emissionResidueMap["CC-"]= 34;
	emissionResidueMap["CTA"]= 35;
	emissionResidueMap["CTC"]= 36;
	emissionResidueMap["CTT"]= 37;
	emissionResidueMap["CTG"]= 38;
	emissionResidueMap["CT-"]= 39;
	emissionResidueMap["CGA"]= 40;
	emissionResidueMap["CGC"]= 41;
	emissionResidueMap["CGT"]= 42;
	emissionResidueMap["CGG"]= 43;
	emissionResidueMap["CG-"]= 44;
	emissionResidueMap["C-A"]= 45;
	emissionResidueMap["C-C"]= 46;
	emissionResidueMap["C-T"]= 47;
	emissionResidueMap["C-G"]= 48;
	emissionResidueMap["C--"]= 49;
	emissionResidueMap["TAA"]= 50;
	emissionResidueMap["TAC"]= 51;
	emissionResidueMap["TAT"]= 52;
	emissionResidueMap["TAG"]= 53;
	emissionResidueMap["TA-"]= 54;
	emissionResidueMap["TCA"]= 55;
	emissionResidueMap["TCC"]= 56;
	emissionResidueMap["TCT"]= 57;
	emissionResidueMap["TCG"]= 58;
	emissionResidueMap["TC-"]= 59;
	emissionResidueMap["TTA"]= 60;
	emissionResidueMap["TTC"]= 61;
	emissionResidueMap["TTT"]= 62;
	emissionResidueMap["TTG"]= 63;
	emissionResidueMap["TT-"]= 64;
	emissionResidueMap["TGA"]= 65;
	emissionResidueMap["TGC"]= 66;
	emissionResidueMap["TGT"]= 67;
	emissionResidueMap["TGG"]= 68;
	emissionResidueMap["TG-"]= 69;
	emissionResidueMap["T-A"]= 70;
	emissionResidueMap["T-C"]= 71;
	emissionResidueMap["T-T"]= 72;
	emissionResidueMap["T-G"]= 73;
	emissionResidueMap["T--"]= 74;
	emissionResidueMap["GAA"]= 75;
	emissionResidueMap["GAC"]= 76;
	emissionResidueMap["GAT"]= 77;
	emissionResidueMap["GAG"]= 78;
	emissionResidueMap["GA-"]= 79;
	emissionResidueMap["GCA"]= 80;
	emissionResidueMap["GCC"]= 81;
	emissionResidueMap["GCT"]= 82;
	emissionResidueMap["GCG"]= 83;
	emissionResidueMap["GC-"]= 84;
	emissionResidueMap["GTA"]= 85;
	emissionResidueMap["GTC"]= 86;
	emissionResidueMap["GTT"]= 87;
	emissionResidueMap["GTG"]= 88;
	emissionResidueMap["GT-"]= 89;
	emissionResidueMap["GGA"]= 90;
	emissionResidueMap["GGC"]= 91;
	emissionResidueMap["GGT"]= 92;
	emissionResidueMap["GGG"]= 93;
	emissionResidueMap["GG-"]= 94;
	emissionResidueMap["G-A"]= 95;
	emissionResidueMap["G-C"]= 96;
	emissionResidueMap["G-T"]= 97;
	emissionResidueMap["G-G"]= 98;
	emissionResidueMap["G--"]= 99;
}

// int getIndex(char residue)
//  Purpose: 
//	  Returns the index in the emission probabilities for the residue
int HMMProbabilities::getEmissionResidueIndex(string residue) {
	return emissionResidueMap.at(residue);
}

void HMMProbabilities::populateEmissionProbabilities(int state, string file) {
	ifstream inputFile(file);
	int totalCount = 0;
	string line;

	// Set the probabilty to the count from the file (we will divide by the total
	//  count in the next step to get the actual probability)
	while(getline(inputFile, line)) {
		vector<string> tokens;
		StringUtilities::split(line, '\t', tokens);
		int count = atoi(tokens.back().c_str());
		setEmissionProbability(state, tokens.front(), count);
		totalCount += count;
	}

	//  Done reading in file 
	inputFile.close();

	// Divide the count set in the first step by the total count to get the probability
	for (pair<string, int> mapPair : emissionResidueMap) {
		string& residue = mapPair.first;
		long double probability =
			emissionProbability(state, residue) / (double) totalCount;
		setEmissionProbability(state, residue, probability);
	}

}

