/*
 * HMMNode.h
 *
 *	This is the header file for the HMMNode object. The HMMNode represents a node
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

#ifndef HMMNODE_H
#define HMMNODE_H
#include <vector>
#include <map>
#include <string>
using namespace std;

class HMMTransition;
class HiddenMarkovModel;

class HMMNode
{
public:
	// Constuctors
	// ==============================================
	HMMNode();
	HMMNode(int anId, int aState, string aResidue, HiddenMarkovModel* aModel);

	// Destructor
	// =============================================
	~HMMNode();

	// Public Attributes
	// =============================================
	int id;
	int state;
	string residue;
	vector<HMMTransition*> inTransitions;
	vector<HMMTransition*> outTransitions;
	double highestWeight;
	HMMNode* highestWeightPreviousNode;
	long double logForwardProbability;
	long double logBackwardProbability;
	long double logConditionalProbability;
	HiddenMarkovModel* model;

	// Public Methods
	// =============================================

	// addInTransition(HMMTransition* aTranstion)
	//  Purpose: 
	//		Add aTransition to the inTranstions collection	
	//  Postconditions:
	//		inTransitions - contains aTransition
	void addInTransition(HMMTransition* aTranstion);

	// addOutTransition(HMMTransition* aTranstion)
	//  Purpose: 
	//		Add aTransition to the outTranstions collection	
	//  Postconditions:
	//		outTransitions - contains aTransition
	void addOutTransition(HMMTransition* aTranstion);

	// 	double logEmissionProbability()
	//  Purpose: 
	//		Returns the log of the emission probability for the state and
	//		residue on this node
	long double logEmissionProbability();

};

#endif // HMMNODE_H
