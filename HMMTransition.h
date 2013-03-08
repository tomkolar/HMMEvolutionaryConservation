/*
 * HMMTransition.h
 *
 *	This is the header file for the HMMTransition object. HMMTransition
 *  represents a transtion for a hidden markov model.  It simply consists
 *  of a start node and end node.  State information from each node can be
 *  used to determine what type of transition this is (same state to same
 *  state? state1 to state2, state2 to state1, etc.).
 *
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */

#ifndef HMMTRANSITION_H
#define HMMTRANSITION_H
#include "HMMNode.h"
using namespace std;

class HMMTransition
{
public:
	// Constuctors
	// ==============================================
	HMMTransition();
	HMMTransition(HMMNode* aStartNode, HMMNode* anEndNode, HiddenMarkovModel* aModel);

	// Destructor
	// =============================================
	~HMMTransition();

	// Public Attributes
	// =============================================
	HMMNode* startNode;
	HMMNode* endNode;
	long double logConditionalProbability;
	HiddenMarkovModel* model;

	// Public Methods
	// =============================================

	// double logProbability()
	//  Purpose: 
	//		Returns the log of the transition probability for transitioning
	//		from the startNode to the endNode
	long double logProbability();
};

#endif // HMMTRANSITION_H
