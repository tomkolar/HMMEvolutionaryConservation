#include "HMMPosition.h"
#include "HMMTransition.h"
#include "MathUtilities.h"
#include <stddef.h>
#include <limits>
#include <math.h>

// Constuctors
// ==============================================
HMMPosition::HMMPosition() {
	// 
	id = 0;
	HMMNode* startNode = new HMMNode();
	nodes.push_back(startNode);
}

HMMPosition::HMMPosition(int anId, string residue, int numStates, HiddenMarkovModel* model) {
	id = anId;
	for (int state = 1; state < numStates; state++) {
		HMMNode* node = new HMMNode(anId, state, residue, model);
		nodes.push_back(node);
	}
}

// Destructor
// =============================================
HMMPosition::~HMMPosition() {
}

// Public Methods
// =============================================

// HMMNode* highestScoringNode()
//  Purpose: 
//		Returns the node with the highest weight from the
//		collection of nodes.
HMMNode* HMMPosition::highestScoringNode() {
	HMMNode* highestScorer = NULL;
	for (HMMNode* node : nodes) {
		if (highestScorer == NULL)
			highestScorer = node;
		else {
			if (node->highestWeight > highestScorer->highestWeight) 
				highestScorer = node;
		}
	}

	return highestScorer;
}

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
void HMMPosition::calculateLogForwardProbability() {

	// Skip start node
	if (id == 0)
		return;

	// Calculation for first node only
	if (id == 1) {
		for (HMMNode* node : nodes) {
			node->logForwardProbability = 
				MathUtilities::elnprod(
					node->inTransitions[0]->logProbability(),    // Initiation prob
					node->logEmissionProbability()				 // Emission prob
				);
		}
		return;
	}

	// Calculation for all other nodes
	for (HMMNode* node : nodes) {
		long double logAlpha = std::numeric_limits<double>::quiet_NaN();
		for (HMMTransition* inTransition : node->inTransitions) {
			logAlpha = 
				MathUtilities::elnsum(
					logAlpha,
					MathUtilities::elnprod(
						inTransition->startNode->logForwardProbability,		// prev prob
						inTransition->logProbability()						// transition prob
					)
				);
		}
		node->logForwardProbability = 
			MathUtilities::elnprod(
				logAlpha,
				node->logEmissionProbability()								 // emission prob
			);
	}
}

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
void HMMPosition::calculateLogBackwardProbability() {

	for (HMMNode* node : nodes) {
		long double logBeta = std::numeric_limits<double>::quiet_NaN();
		for (HMMTransition* outTransition : node->outTransitions) {
			logBeta =
				MathUtilities::elnsum(
					logBeta,
					MathUtilities::elnprod(
						outTransition->logProbability(),						// transition prob
						MathUtilities::elnprod(
							outTransition->endNode->logEmissionProbability(),	// emission prob
							outTransition->endNode->logBackwardProbability		// prev prob
						)
					)
				);
		}
		node->logBackwardProbability = logBeta;	
	}
}

// calculateLogConditionalProbability()
//  Purpose: 
//		Calculate and store the log conditional probability for all the nodes.
//		This is the probabilty of being in a state at a particular at this position
//		given the model.
//
//  Postconditions:
//		nodes.logConditionalProbabilty will be set for all nodes
void HMMPosition::calculateNodeLogConditionalProbabilities() {

	// Calculated normailzer for forwardProp*backwardProb for this position
	//   (Sum up forwardProb*backwardProb for all nodes) 
	long double normalizer = std::numeric_limits<double>::quiet_NaN();
	for (HMMNode* node : nodes) {
		// Skip start node
		if (node->id == 0) {
			continue;
		}

		node->logConditionalProbability =
			MathUtilities::elnprod(
				node->logForwardProbability,
				node->logBackwardProbability
			);
		normalizer =
			MathUtilities::elnsum(normalizer, node->logConditionalProbability);
	}

	// Calculate the condtional probability for each node at this position
	//	(forwardProb*backwardProp/normalizer)
	for (HMMNode* node : nodes) {
		// Skip start node
		if (node->id == 0) {
			continue;
		}

		node->logConditionalProbability =
			MathUtilities::elnprod(
				node->logConditionalProbability,
				-normalizer
			);
	}

}

void HMMPosition::calculateTransitionLogConditionalProbabilities() {

	// Calculate normailzer and non-normalized log conditional prob
	// for each transition
	long double normalizer = std::numeric_limits<double>::quiet_NaN();
	for (HMMNode* node : nodes) {
		for (HMMTransition* outTransition : node->outTransitions) {
		
			outTransition->logConditionalProbability =
				MathUtilities::elnprod(
					node->logForwardProbability,							// forward prob
					MathUtilities::elnprod(
						outTransition->logProbability(),						// transition prob
						MathUtilities::elnprod(
							outTransition->endNode->logEmissionProbability(),	// next emission prob
							outTransition->endNode->logBackwardProbability   // next backward prob
						)
					)
				);
			normalizer =
				MathUtilities::elnsum(normalizer, outTransition->logConditionalProbability);
		}
	}

	// Normalize the cacluated values
	for (HMMNode* node : nodes) {
		for (HMMTransition* outTransition : node->outTransitions) {
		
			outTransition->logConditionalProbability =
				MathUtilities::elnprod(
					outTransition->logConditionalProbability,
					-normalizer
				);
		}
	}
}

double HMMPosition::logLikelihood() {

	double logLikelihood = std::numeric_limits<double>::quiet_NaN();

	for (HMMNode* node : nodes) {
		logLikelihood = MathUtilities::elnsum(logLikelihood,node->logForwardProbability);
	}

	return logLikelihood / log(2);
} 

double HMMPosition::logLikelihoodBackward() {

	double logLikelihood = std::numeric_limits<double>::quiet_NaN();

	for (HMMNode* node : nodes) {
		logLikelihood = MathUtilities::elnsum(logLikelihood,node->logBackwardProbability);
	}

	return logLikelihood / log(2);
} 
