/*
 * HMMNode.cpp
 *
 *	This is the cpp file for the HMMNode object. The HMMNode represents a node
 *  in a hidden markov model.  Mostly it is a place holder for information
 *  pertaining to a particular state at a position in the HMM.
 *
 *  Important Attributes:
 *		id - position in the HMM
 *		state - the underlying state this node represents in the HMM
 *		residue - the residue in this position in the sequence used to build the HMM
 *		inTransitions - HMMTransitions that are coming into this node from the previous
 *						position in the HMM
 *		outTransitions - HMMTransitions that are leaving this node to the next position
 *						 in the HMM
 *		highestWeight - the highest weight determined by the viterbi path
 *		higheestWeightPreviousNode - the previous node in the HMM that gave the highest
 *									 weight viterbi path
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */
#include "HMMNode.h"
#include "HiddenMarkovModel.h"
#include "HMMProbabilities.h"

// Constuctors
// ==============================================
HMMNode::HMMNode(){
	// Initialize as start node
	id = 0;
	state = 0;
	residue = "";
	highestWeight = 0; // set to zero only on start node
	highestWeightPreviousNode = NULL;
	logForwardProbability = 0;
	logBackwardProbability = 0;
	logConditionalProbability = 0;
}

HMMNode::HMMNode(int anId, int aState, string aResidue, HiddenMarkovModel* aModel) {
	id = anId;
	state = aState;
	residue = aResidue;
	model = aModel;
	highestWeightPreviousNode = NULL;
}

// Destructor
// =============================================
HMMNode::~HMMNode()
{
}

// Public Methods
// =============================================

// addInTransition(HMMTransition* aTranstion)
//  Purpose: 
//		Add aTransition to the inTranstions collection	
//  Postconditions:
//		inTransitions - contains aTransition
void HMMNode::addInTransition(HMMTransition* aTransition) {
	inTransitions.push_back(aTransition);
}

// addOutTransition(HMMTransition* aTranstion)
//  Purpose: 
//		Add aTransition to the outTranstions collection	
//  Postconditions:
//		outTransitions - contains aTransition
void HMMNode::addOutTransition(HMMTransition* aTransition) {
	outTransitions.push_back(aTransition);
}

// 	double logEmissionProbability()
//  Purpose: 
//		Returns the log of the emission probability for the state and
//		residue on this node
long double HMMNode::logEmissionProbability() {
	return model->probabilities->logEmissionProbability(state, residue);
}

