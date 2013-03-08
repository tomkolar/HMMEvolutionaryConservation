/*
 * driver.cpp
 *
 *	This is the driver file for creating an a hidden markov model from
 *  a fastafile.  Viterbi training is then exectued to try and locate
 *  GC rich portions of the sequence.
 *
 *	Typical use:
 *		hmm fastaFile numIterations
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */
#include "MultipleAlignmentFile.h"
#include "HiddenMarkovModel.h"
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

int main( int argc, char *argv[] ) {
/*
	// Check that file name and iterations were entered as arguments
	if (argc < 5) {
			cout << "Invalid # of arguments\n";
			cout << "usage: hmm multipleAlignmentFile iterations nuetralCountsFile conservedCountsFile \n";
			return -1;
	}

	cout << "Starting\n";

	// Get Parameters
	string multiAlignFileName = argv[1];
	int iterations = atoi(argv[2]);
	string neutralCountsFileName = argv[3];
	string conservedCountsFileName = argv[4];
*/
	// Set Parameters
	string multiAlignFileName = "c:/Users/kolart/Documents/Genome540/Assignment8/ENm006_short.aln";
	int iterations = 1;
	string neutralCountsFileName = "c:/Users/kolart/Documents/Genome540/Assignment8/anc_rep_counts.txt";
	string conservedCountsFileName = "c:/Users/kolart/Documents/Genome540/Assignment8/codon1_2_counts.txt";

	// Create the fasta file object
	MultipleAlignmentFile* multiAlignFile =
		new MultipleAlignmentFile(multiAlignFileName);
	cout << "Multi Align's done\n";

	// Create the Hidden Markov Model
	HiddenMarkovModel hmm(multiAlignFile, neutralCountsFileName, conservedCountsFileName);
	hmm.viterbiTraining(iterations);

	cout << hmm.viterbiResultsString();
	cout << "Fred";
}