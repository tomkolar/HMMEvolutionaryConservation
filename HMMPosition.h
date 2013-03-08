/*
 * HMMPosition.h
 *
 *	This is the header file for the HMMPosition object. The HMMPoistion
 *  represents a position in a hidden markov model.  The postion is
 *  essentially a collection of HMMNodes.  One node for each state
 *  in the HMM.
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */

#ifndef HMMPOSITION_H
#define HMMPOSITION_H
#include "HMMNode.h"
using namespace std;

class HiddenMarkovModel;

class HMMPosition
{
public:
	// Constuctors
	// ==============================================
	HMMPosition();
	HMMPosition(int anId, string residue, int numStates, HiddenMarkovModel* model);

	// Destructor
	// =============================================
	~HMMPosition();

	// Public Attributes
	// =============================================
	vector<HMMNode*> nodes;
	int id;

	// Public Methods
	// =============================================

	// HMMNode* highestScoringNode()
	//  Purpose: 
	//		Returns the node with the highest weight from the
	//		collection of nodes.
	HMMNode* highestScoringNode();

	// calculateLogForwardProbabilty()
	//  Purpose: 
	//		Calculate and store the log forward probabilty for the forward-backward
	//		(Baum-Welch) algorithm.
	//
	//		The forward probability determined seperately for each node.  The forward
	//		probability for a particular node is determined by itetarting through each
	//		of the and summing the calculated forward probability from each.  The calculation
	//		for each node is the same as the formula for highest weight path:
	//			  previous nodes weight
	//			+ transmission probability
	//			+ emission Probability
	//
	//  Postconditions:
	//		logForwardProb - set to calculated log probability
	void calculateLogForwardProbability();

	// calculateLogBackwardProbabilty()
	//  Purpose: 
	//		Calculate and store the log backward probabilty for the forward-backward
	//		(Baum-Welch) algorithm.
	//
	//		The backward probability determined seperately for each node.  The forward
	//		probability for a particular node is determined by itetarting through each
	//		of the and summing the calculated forward probability from each.  The calculation
	//		for each node is the same as the formula for highest weight path:
	//			  previous nodes weight
	//			+ transmission probability
	//			+ emission Probability
	//
	//  Postconditions:
	//		logBackwwardProb - set to calculated log probability
	void calculateLogBackwardProbability();

	// calculateLogConditionalProbability()
	//  Purpose: 
	//		Calculate and store the log conditional probability for all the nodes.
	//		This is the probabilty of being in a state at a particular at this position
	//		given the model.
	//
	//  Postconditions:
	//		nodes.logConditionalProbabilty will be set for all nodes
	void calculateNodeLogConditionalProbabilities();

	void calculateTransitionLogConditionalProbabilities();
	double logLikelihood();
	double logLikelihoodBackward();

};

#endif // HMMPOSITION_H
