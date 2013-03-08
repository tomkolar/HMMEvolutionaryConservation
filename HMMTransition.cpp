/*
 * HMMTransition.cpp
 *
 *	This is the cpp file for the HMMTransition object. HMMTransition
 *  represents a transtion for a hidden markov model.  It simply consists
 *  of a start node and end node.  State information from each node can be
 *  used to determine what type of transition this is (same state to same
 *  state? state1 to state2, state2 to state1, etc.).
 *
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */
#include "HMMTransition.h"
#include "HiddenMarkovModel.h"

// Constuctors
// ==============================================
HMMTransition::HMMTransition() {
}

HMMTransition::HMMTransition(HMMNode* aStartNode, HMMNode* anEndNode, HiddenMarkovModel* aModel) {
	startNode = aStartNode;
	endNode = anEndNode;
	model = aModel;
}

// Destructor
// =============================================
HMMTransition::~HMMTransition() {
}

// Public Methods
// =============================================

// double logProbability()
//  Purpose: 
//		Returns the log of the transition probability for transitioning
//		from the startNode to the endNode
long double HMMTransition::logProbability() {
	if (startNode->state == 0)
		return model->probabilities->logInitiationProbability(endNode->state);
	else
		return model->probabilities->logTransitionProbability(startNode->state, endNode->state);
}
