/*
 * HiddenMarkovModel.h
 *
 *	This is the header file for the HiddenMarkovModel object. 
 *  HiddenMarkovModel is the main object representing a hidden
 *  markov model.  At the moment it is hard coded to support
 *  two states.  However, it can easily be modified to support
 *  more states (see details at bottom of this header comment).
 *
 *	The model attribute holds the generated hidden markov model.
 *  Essentially the model is a vector of HMMPosition objects. A
 *  HMMPosition object is basically a collection of HMMNodes (one for
 *  each state in the HMM).  The nodes know there position, state,
 *  and residue and contain collections of HHMTransitions that come
 *  in to the node and out from the node. Emission probabilites can be 
 *  determined from the state/residue information.  Transition probabilites
 *  are determined by the start/stop node states for an HMMTransition
 *  object.
 *
 *	The probabilitites attribute holds the inititation, emission and
 *  transition probabilties that are used when finding a path through
 *  the model.  HMMNodes will access the probabilites to determine their
 *  emission probabilities.  HMMTransitions will access the probabilties
 *  to determine their initiation/transition probabilties.
 *
 *  Viterbi training is currently the only implmented method for creating
 *  a path.  Typical use would be:
 *
 *		HiddenMarkovModel(aFastaFile)
 *			- instantate the object with the fasta file to build the
 *			  model from
 *
 *		viterbiTraining(numIterations)
 *			- run viterbi training to generate the model and then train
 *			  for the number of iterations specified
 *
 *		viterbiResultsString()
 *			- returns a string with the results from each iteration of the
 *			  viterbi training
 *
 *	Additional Methods of interest:
 *
 *		allScoresResultsString()
 *			- returns a string with the score for each state at each position
 *			  after the last run of  viterbi training
 *
 *		pathStatesResultsString()
 *			- returns a string of the state for every position in the 
 *			  viterbi path
 *
 *  Reconfiguration for more states:
 *	  The ultimate fix would be to take in the number of states
 *    as a parameter.  There is already a numStates variable that
 *    exists and is hardcoded to 2, but this could be set to the
 *    passed in parameter.  The main issue is with the HMMProbabilites
 *    object as it contains several arrays that have sizes the same
 *    size as the number of states.  These would need to be initialized
 *    to have the correct size, or perhaps changed to vectors (as it
 *    is easier to initialize vectors to a particular size).
 *
 *  Created on: 2-13-13
 *      Author: tomkolar
 */

#ifndef HIDDENMARKOVMODEL_H
#define HIDDENMARKOVMODEL_H
#include "MultipleAlignmentFile.h"
#include "HMMPosition.h"
#include "HMMProbabilities.h"
#include "HMMViterbiResults.h"
#include <vector>
#include <map>
using namespace std;

class HiddenMarkovModel
{
public:
	// Constuctors
	// ==============================================
	HiddenMarkovModel();
	HiddenMarkovModel(MultipleAlignmentFile* aMultiAlginFile,
		string neutralCountsFileName, string conservedCountsFileName);

	// Destructor
	// =============================================
	~HiddenMarkovModel();

	// Public Attributes
	// =============================================
	HMMProbabilities* probabilities;
	vector<HMMViterbiResults*> viterbiResults;

	// Public Methods
	// =============================================

	// viterbiTraining(int numIterations)
	//  Purpose: 
	//		Perform viterbi training for the number of iterations specified.
	//
	//		Each iteration consists of the following steps
	//			1. Build Hidden Markov Model and calculate viterbi weight
	//			   for each node
	//			2. Walk the viterbi path backwards and gather the 
	//			   results (HMMViterbiResults)
	//			3. Reset the probabilites to the ones calculated from the'
	//			   Viterbi results
	//
	//  Postconditions:
	//		viterbiResults - contains results from each training iteration
	//		probabilities - modified at end of each training iteration to 
	//						reflect calculated probabilities from the viterbi
	//						results
	void viterbiTraining(int numIterations);

	// baumWelchTraining()
	//  Purpose: 
	//		Use the Baum-Welch (forward-backward) algorithm to estimate
	//		the paramters for the model
	//
	//		Each iteration consists of the following steps
	//			1. Build Hidden Markov Model and calculate forward-backward weight
	//			   for each node
	void baumWelchTraining();

	// string allScoresResultsString()
	//  Purpose:
	//		Returns a string representing the score (weight) from each node
	//		in each position of the viterbi path.
	//
	//		format:
	//			Position: <positionId>
	//			  Node: (<node1State>,<node1Weight>)
	//			  Node: (<node2State>,<node2Weight>)
	//			  ...
	//  Preconditions:
	//		viterbiTraining has been run
	string allScoresResultsString();

	// string pathStatesResultsString()
	//  Purpose:
	//		Returns a string representing the state for each position in the
	//		viterbi path.
	//
	//		format:
	//			<position1State><position2State> ... <positionNstate>
	//  Preconditions:
	//		viterbiTraining has been run
	string pathStatesResultsString();

	// string viterbiResultsString()
	//  Purpose:
	//		Returns a string representing the results for each iteration in
	//		the viterbi training. See HMMViterbiResults.resultsWithoutSegments()
	//		and HMMViterbiResults.allResults() for details on the format of
	//		the viterbi results for a particular iteration.
	//
	//		format:
	//			<viterbiResultsIteration1.resultsWithouSegments()>
	//			<viterbiResultsIteration2.resultsWithouSegments()>
	//			...
	//			<viterbiResultsIterationLast.allResults()>
	//  Preconditions:
	//		viterbiTraining has been run
	string viterbiResultsString();

private:

	// Private Attributes
	// =============================================
	static const int numStates;
	MultipleAlignmentFile* multiAlignFile;
	vector<HMMPosition*> model;
	bool modelBuilt;

	// Private Methods
	// =============================================

	// buildAndCalculateModel(bool calculateForward)
	//  Purpose: 
	//		Build the hidden markov model (if not already built) and calculate
	//		the viterbi weight or the forward probability for each node.
	//
	//		If the model has already been built, then this method will simply
	//		(re)calcuate the viterbi weight or forward probability for each node. 
	//		The expectation is that the probabilities have been set to new values
	//		and we are recalculating the weights using these new probabilities.
	//
	//		Model Building Steps:
	//			1. Create Start Node
	//			2. Iterate through the sequence from the fasta file and do the
	//			   following:
	//				a. create a HMMPosition object with one node for each state
	//				b. create HMMTransitions for this position
	//				c. calculate the viterbi weight or forwared probability for each
	//				   node in the position
	//				d. add the position to the model attribute
	//
	//  Postconditions:
	//		model - contains HMMPosition objects for every position in the
	//				sequence from the fastaFile
	void buildAndCalculateModel(bool calculateForward);

	// createTransitionsFor(HMMPosition* currentPosition, HMMPosition* previousPosition)
	//  Purpose: 
	//		Create HMMTransitions from the previousPosition to the currentPosition.
	//		Transition objects wlill be created for each node in the previousPostion
	//		to each node in the currentPosition
	//  Postconditions:
	//		currentPosition.inTransitions
	//			- contains transitions from each node in the previous position
	//		previousPosition.outTransitions
	//			- contains transitions to each node in the current position
	void createTransitionsFor(HMMPosition* currentPosition, HMMPosition* previousPosition);

	// calculateHighestWeightPath(HMMPosition* currentPosition)
	//  Purpose: 
	//		Calculate and store the highest weight using the viterbi algorithm
	//		(the incoming path with the highest weight).
	//
	//		Highest Weight is determined seperately for each node.  Highest wieght
	//		for a particular node is determined by itetarting through each of the 
	//		incoming transitions taking the one with the highest weight.  The path
	//		weight uses the formula:
	//			  previous nodes weight
	//			+ transmission probability
	//			+ emission Probability
	//
	//  Postconditions:
	//		currentPosition.highestWeight - set to highest calculated weight
	//		currentPosition.highestWeightPreviousNode
	//			- set to the previous node that generated the highest calculated weight
	void calculateHighestWeightPath(HMMPosition* currentPosition);

	// calculateLogBackwaredPorbabilities()
	//  Purpose: 
	//		Calculate and store the log backwared probability for all the nodes
	//		in the model.
	//
	//  Postconditions:
	//		model.nodes.logBackwardProbabilty will be set for all nodes
	void calculateLogBackwardProbabilities();

	// calculateLogConditionalProbabilities()
	//  Purpose: 
	//		Calculate and store the log conditional probability for all the nodes
	//		in the model.  This is the probabilty of being in a state at a particular
	//		position given the model.
	//
	//  Postconditions:
	//		model.nodes.logConditionalProbabilty will be set for all nodes
	void calculateLogConditionalProbabilities();

	// HMMViterbiResults* gatherViterbiResults(int iteration);
	//  Purpose: 
	//		Creates, populates, and return a HMMViterbiResults object containing
	//		the results for the most recent iteration in the viterbi training.
	//
	//		Results are gathered by walking the viterbi path backward.  Results
	//		gathered include the following:
	//			state counts - how many times a state occurs in the path
	//			segment counts - how many segments (i.e., continuos occurencee of
	//							 one state) of each state occurs
	//			segments - collection of start and stop values for each segment
	//			transition counts - counts for how many time states transition (both
	//								from one state to another and from one state to 
	//							    the stame state)
	//
	//		After the results are gathered, the probabitlites will be recalculated
	//		using the gathered information.
	//
	//  Postconditions:
	//		currentPosition.highestWeight - set to highest calculated weight
	//		currentPosition.highestWeightPreviousNode
	//			- set to the previous node that generated the highest calculated weight
	HMMViterbiResults* gatherViterbiResults(int iteration);

	void calculateBaumWelchEmissionProbabilities();
	void calculateBaumWelchTransitionProbabilities();
	void calculateBaumWelchInitiationProbabilities();
	double calculateLogLikelihood();
	string baumWelchResultsString(int iterations, double logLikelihood);


};

#endif // HIDDENMARKOVMODEL_H
