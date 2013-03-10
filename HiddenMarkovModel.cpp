/*
 * HiddenMarkovModel.cpp
 *
 *	This is the cpp file for the HiddenMarkovModel object. 
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

#include "HiddenMarkovModel.h"
#include "HMMTransition.h"
#include "HMMProbabilities.h"
#include "MathUtilities.h"
#include <sstream>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <limits>
#include <stdexcept> 

// const variable initialization
// ==============================================
const int HiddenMarkovModel::numStates = 3;

// Constuctors
// ==============================================
HiddenMarkovModel::HiddenMarkovModel() {
}

HiddenMarkovModel::HiddenMarkovModel(MultipleAlignmentFile* aMultiAlignFile,
		string neutralCountsFileName, string conservedCountsFileName) {
	multiAlignFile = aMultiAlignFile;
	modelBuilt = false;
	probabilities = HMMProbabilities::initialProbabilities(neutralCountsFileName, conservedCountsFileName);
	
	// Print out the initial probabilities
	cout << probabilities->probabilitiesResultsString();
}

// Destructor
// =============================================
HiddenMarkovModel::~HiddenMarkovModel() {
}

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
void HiddenMarkovModel::viterbiTraining(int numIterations) {

	for (int iteration = 1; iteration <= numIterations; iteration++) {
		// Build the model and calculate the weights
		buildAndCalculateModel(false);
		cout << "Model Built.\n";

		// Gather the viterbi reuslts
		HMMViterbiResults* aViterbiResults = gatherViterbiResults(iteration);
		viterbiResults.push_back(aViterbiResults);
		cout << "Viterbi Results Gathered.\n";

		// Reset the probabilities to the viterbi calculated ones for the next
		// iteration
		probabilities = aViterbiResults->probabilities;
	}
}

// baumWelchTraining()
//  Purpose: 
//		Use the Baum-Welch (forward-backward) algorithm to estimate
//		the paramters for the model
//
//		Each iteration consists of the following steps
//			1. Build Hidden Markov Model and calculate forward-backward weight
//			   for each node
void HiddenMarkovModel::baumWelchTraining() {
	bool trainingDone = false;
	int iterationCounter = 0;
	double previousLogLikelihood = 0;
	while (!trainingDone) {
		// Build the model and calculate the forward/backward probabilites
		buildAndCalculateModel(true);
		calculateLogBackwardProbabilities();
		calculateLogConditionalProbabilities();

		// Calculate the new transition/emission probabilties
		calculateBaumWelchEmissionProbabilities();
		calculateBaumWelchInitiationProbabilities();
		calculateBaumWelchTransitionProbabilities();

		// Calcualte likelihood and check if done
		double currentLogLikelihood = model.back()->logLikelihood();
		if (abs(previousLogLikelihood - currentLogLikelihood) < 0.1)
			trainingDone = true;

		// Set values for next iteration
		previousLogLikelihood = currentLogLikelihood;
		iterationCounter++;
		cout
			<< "Iteration: " << iterationCounter 
			<< "  Likelihood: " << currentLogLikelihood
			<< "\n";
	}

	cout << baumWelchResultsString(iterationCounter, previousLogLikelihood);
}

string HiddenMarkovModel::baumWelchResultsString(int iterations, double logLikelihood) {
	stringstream ss;

	// EM Result header
	ss << "    <result type=\"EM_result\">\n";

	// Iterations
	ss
		<< "      <result type=\"iterations\">"
		<< iterations
		<< "</result>\n";

	// Log Likelihood
	ss
		<< "      <result type=\"log_likelihood\">"
		<< logLikelihood
		<< "</result>\n";

	// Probabilities
	ss << probabilities->probabilitiesResultsString();

	// EM Result footer
	ss << "    </result>\n";

	return ss.str();
}

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
string HiddenMarkovModel::allScoresResultsString() {
	stringstream ss;

	for (HMMPosition* aPosition : model) {
		ss << "Position: " << aPosition->id << "\n";
		for (HMMNode* node : aPosition->nodes) {
			ss << "  Node: ("
			   <<  node->state
			   << ", "
			   << node->highestWeight
			   << ")\n";
		}
	}

	return ss.str();
}

// string pathStatesResultsString()
//  Purpose:
//		Returns a string representing the state for each position in the
//		viterbi path.
//
//		format:
//			<position1State><position2State> ... <positionNstate>
//  Preconditions:
//		viterbiTraining has been run
string HiddenMarkovModel::pathStatesResultsString() {
	stringstream ss;

	HMMPosition* lastPosition = model.back();
	HMMNode* aNode = lastPosition->highestScoringNode();
	while (aNode->state != 0) {
		ss << aNode->state;
		aNode = aNode->highestWeightPreviousNode;
	}

	string reversePath = ss.str();
	return string(reversePath.rbegin(), reversePath.rend());
}

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
string HiddenMarkovModel::viterbiResultsString() {
	stringstream ss;

	int numResults = viterbiResults.size();
	for (int i = 0; i < numResults; i++) {
		if (i < numResults -1)
			ss << viterbiResults[i]->resultsWithoutSegments();
		else
			ss << viterbiResults[i]->allResults();
	}
	
	return ss.str();
}

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
void HiddenMarkovModel::buildAndCalculateModel(bool calculateForward) {

	if (modelBuilt) {
		// Model already built so just recalculate the weights using
		// the newly set probabilities
		for (HMMPosition* aPosition : model) {
			if (calculateForward)
				aPosition->calculateLogForwardProbability();
			else
				calculateHighestWeightPath(aPosition);
		}
	}
	else {
		vector<string>& sequence = multiAlignFile->getSequence();

		// Create Start Position
		HMMPosition* startPosition = new HMMPosition();
		model.push_back(startPosition);

		// Iterate through the sequence and create model on the fly
		HMMPosition* previousPosition = startPosition;
		int seqLength = sequence.size();
		for (int seqPos = 0; seqPos <= seqLength - 1; seqPos++) {
			// Create a Position object with one node for each state
			HMMPosition* aPosition = new HMMPosition(seqPos, sequence[seqPos], numStates, this);

			// Create the incoming transitions for the curent position
			createTransitionsFor(aPosition, previousPosition);

			// Calculate forward probability or highest weight path
			if (calculateForward)
				aPosition->calculateLogForwardProbability();
			else
				calculateHighestWeightPath(aPosition);

			// Add position to the model
			model.push_back(aPosition);

			// Set previous positon to current position so that next iteration cancreate
			// transitions correctly.
			previousPosition = aPosition;
		}

		modelBuilt = true;
	}
}

// calculateLogBackwaredPorbabilities()
//  Purpose: 
//		Calculate and store the log backwared probability for all the nodes
//		in the model.
//
//  Postconditions:
//		model.nodes.logBackwardProbabilty will be set for all nodes
void HiddenMarkovModel::calculateLogBackwardProbabilities() {
	int numPositions = model.size();

	// Set the backward probabilty on the last position equal to 0
	HMMPosition* lastPosition = model[numPositions - 1];
	for (HMMNode* node : lastPosition->nodes) {
		node->logBackwardProbability = 0.0;
	}

	// Walk the positions backward and calculate the probabilites
	for (int positionId = numPositions - 2; positionId >= 1; positionId--) {
		HMMPosition* position = model[positionId];
		position->calculateLogBackwardProbability();
	}
}

// calculateLogConditionalProbabilities()
//  Purpose: 
//		Calculate and store the log conditional probability for all the nodes
//		in the model.  This is the probabilty of being in a state at a particular
//		position given the model.
//
//  Postconditions:
//		model.nodes.logConditionalProbabilty will be set for all nodes
void HiddenMarkovModel::calculateLogConditionalProbabilities() {

	int numPositions = model.size();
	for (HMMPosition* position : model) {
		// Skip start node
		if (position->id == 0)
			continue;

		// Node Probabilities (Gamma)
		position->calculateNodeLogConditionalProbabilities();

		// Transition Probabilites (Epsilon)
		if (position->id != numPositions)  // calculate for all but last node
			position->calculateTransitionLogConditionalProbabilities();
	}
}

void HiddenMarkovModel::calculateBaumWelchEmissionProbabilities() {
/*
	// Create and initialize vectors to track the numerator and denominator
	// calculating the probabilities
	vector<map<char, long double>> numerators;
	vector<long double> denominators;
	for (int i = 0; i < numStates; i++) {

		map<char, long double> numeratorMap;
		numeratorMap['A'] = std::numeric_limits<double>::quiet_NaN();
		numeratorMap['C'] = std::numeric_limits<double>::quiet_NaN();
		numeratorMap['G'] = std::numeric_limits<double>::quiet_NaN();
		numeratorMap['T'] = std::numeric_limits<double>::quiet_NaN();
		numerators.push_back(numeratorMap);
		denominators.push_back(std::numeric_limits<double>::quiet_NaN());
	}

	// Iterate through the positions and populate the numerators and 
	// denominators
	for (HMMPosition* position : model) {
		// Do not calculate for start nodes
		if (position->id == 0)
			continue;

		for (HMMNode* node : position->nodes) {
			// Set numerator for the nodes state and residue
			numerators[node->state].at(node->residue) = 
				MathUtilities::elnsum(
					numerators[node->state].at(node->residue),
					node->logConditionalProbability
				);

			// Set denominator for the nodes state
			denominators[node->state] =
				MathUtilities::elnsum(
					denominators[node->state],
					node->logConditionalProbability
				);
		}
	}

	// Reset the emission probabilities from the calculated numerators
	// and denominators
	for (int state = 0; state < numStates; state++) {
		for (pair<const char, long double> mapPair : numerators[0]) {
			char residue = mapPair.first;
			long double newEmissionProbability = 
				MathUtilities::eexp(
					MathUtilities::elnprod(
						numerators[state].at(residue),
						-denominators[state]
					)
				);
			probabilities->setEmissionProbability(state, residue, newEmissionProbability);
		}
	}
*/
}

void HiddenMarkovModel::calculateBaumWelchInitiationProbabilities(){
	HMMPosition* firstPosition = model[1];
	for (HMMNode* node : firstPosition->nodes) {
		long double newInitiationProbability = MathUtilities::eexp(node->logConditionalProbability);
		probabilities->setInitiationProbability(node->state, newInitiationProbability);
	}
}

void HiddenMarkovModel::calculateBaumWelchTransitionProbabilities() {

	// Create and initialize arrays to track the numerator and denominator
	// calculating the probabilities
	long double numerators[numStates][numStates];
	long double denominators[numStates][numStates];
	for (int i = 0; i < numStates; i++) {
		for (int j = 0; j < numStates; j++) {
			numerators[i][j] = std::numeric_limits<double>::quiet_NaN();
			denominators[i][j] = std::numeric_limits<double>::quiet_NaN();
		}
	}

	// Iterate through the positions and populate the numerators and 
	// denominators
	int numPositions = model.size();
	for (HMMPosition* position : model) {
		// Do not calculate for start nodes or last node
		if (position->id == 0 || position->id == numPositions)
			continue;
		
		for (HMMNode* node : position->nodes) {
			for (HMMTransition* outTransition : node->outTransitions) {
				int startState = outTransition->startNode->state;
				int endState = outTransition->endNode->state;
				numerators[startState][endState] =
					MathUtilities::elnsum(
						numerators[startState][endState],
						outTransition->logConditionalProbability
					);
				denominators[startState][endState] =
					MathUtilities::elnsum(
						denominators[startState][endState],
						node->logConditionalProbability
					);
			}
		}
	}

	// Reset the transition probabilities
	for (int i = 0; i < numStates; i++) {
		for (int j = 0; j < numStates; j++) {
			long double newTransitionProbability =
				MathUtilities::eexp(
					MathUtilities::elnprod(
						numerators[i][j],
						-denominators[i][j]
					)
				);

			probabilities->setTransitionProbability(i,j, newTransitionProbability);
		}
	}
}

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
void HiddenMarkovModel::createTransitionsFor(HMMPosition* currentPosition, HMMPosition* previousPosition) {
	// Create inTransitions for the  current position nodes
	for (HMMNode* currentPositionNode : currentPosition->nodes) {
		// Create one transition for each node from the previous position
		// to this node in the current position
		for (HMMNode* previousPositionNode : previousPosition->nodes) {
			// Create the transition
			HMMTransition* aTransition = new HMMTransition(previousPositionNode, currentPositionNode, this);

			// Add to the transitions collections on nodes
			currentPositionNode->addInTransition(aTransition);
			previousPositionNode->addOutTransition(aTransition);
		}
	}
}

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
void HiddenMarkovModel::calculateHighestWeightPath(HMMPosition* aPosition) {
	// Calculate hwp for each node
	for (HMMNode* positionNode : aPosition->nodes) {
		// initialize highest weight to negative infinity (except start node)
		if (positionNode->id != 0)
			positionNode->highestWeight = -DBL_MAX;

		// Iterater through each of the incoming transitions to find
		// the highest score
		for (HMMTransition* transition : positionNode->inTransitions) {
			// calculate the score
			long double score = 
				MathUtilities::elnprod(
					transition->startNode->highestWeight,  // previous nodes weight
					MathUtilities::elnprod(
						transition->logProbability(),          // transition probability 
						positionNode->logEmissionProbability() // emission probablity
					)
				);

			// Replace the highest weight info if this path has the highest score
			if (!MathUtilities::isNaN(score) && (score > positionNode->highestWeight)) {
				positionNode->highestWeight = score;
				positionNode->highestWeightPreviousNode = transition->startNode;
			}
		}
	}
}
	
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
HMMViterbiResults* HiddenMarkovModel::gatherViterbiResults(int iteration) {
	HMMViterbiResults* results = new HMMViterbiResults(iteration, numStates);

	// Find the end of the highest scoring path
	HMMPosition* lastPosition = model.back();
	HMMNode* aNode = lastPosition->highestScoringNode();

	// Walk the path backward and gather the data
	int previousState = -1; 
	pair<int, int> currentSegment = pair<int,int>(-1,-1);
	int offset = multiAlignFile->getStartPosition();

	while (aNode->id != 0) {
		int currentState = aNode->state;
		
		// Update number of occurrences for a state
		results->stateCounts[currentState]++;

		// Update emission count for state
		results->emissionCounts.at(currentState).at(aNode->residue)++;

		// Update segment info
		if (currentState != previousState) {
			// Set the start of the segment
			currentSegment.first = aNode->id + 1 + offset;

			// Add segment to segments map (except for first time through)
			if (currentSegment.second != -1) 
				results->segments[previousState].push_back(currentSegment);

			// Create new segment
			currentSegment = pair<int,int>(aNode->id + offset, aNode->id + offset);

			// Update number of segments for a state
			results->segmentCounts[currentState]++;
		}

		// Update transition counts
		if (previousState >= 0) {
			results->transitionCounts[currentState][previousState]++;
		}

		// Set up variables for next iteration
		previousState = currentState;
		aNode = aNode->highestWeightPreviousNode;
	}

	// Add the last segment to the collection
	currentSegment.first =  offset;
	results->segments[previousState].push_back(currentSegment);

	// Calculate the probabilities
	results->calculateProbabilities(probabilities);

	// Return the results
	return results;
}
