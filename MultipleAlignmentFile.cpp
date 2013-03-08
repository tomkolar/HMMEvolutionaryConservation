/*
 * MultipleAlignmentFile.h
 *
 *  The MultipleAlignmentFile object is a utility object designed to read in a
 *  Multiple Alignment file and keep the information for the file in memory.
 *
 *  ************* WARNING ***********************************
 *  *  There is no error handling in place for this object! *
 *  *  This means that if the file does not exist or is     *
 *  *  formatted incorrectly, you will get an error and     *
 *  *  will not be able to use this object.                 *
 *  *                                                       *
 *  *  If this code gets moved to a production setting      *
 *  *  appropriate error handling should be implemented!    *
 *  *********************************************************
 *
 * Typical use for the file would be to use the MultipleAlignementFile(fileName)
 * constructor to create the object.  This will automatically open the
 * Multiple Alignment File specified by the fileName, and read its contents
 * storing them in the sequence vector
 *
 *  Created on: 3-6-13
 *      Author: tomkolar
 */

#include "MultipleAlignmentFile.h"
#include "StringUtilities.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

// Constuctors
// ==============================================
MultipleAlignmentFile::MultipleAlignmentFile() {
}

MultipleAlignmentFile::MultipleAlignmentFile(string name) {
	fileName = name;
	populate();
}

// Destructor
// ==============================================
MultipleAlignmentFile::~MultipleAlignmentFile() {
}

// Public Accessors
// =============================================
const int MultipleAlignmentFile::getSequenceLength() {
	return sequence.size();
}

const int MultipleAlignmentFile::getStartPosition() {
	return startPosition;
}

string& MultipleAlignmentFile::getFileName() {
	return fileName;
}

vector<string>& MultipleAlignmentFile::getSequence() {
	return sequence;
}

// Private Methods
// =============================================

// populate()
//  Purpose:
//		Reads in the Multiple Alignment File specified by fileName and
//		populates the object with its contents
//	Preconditions:
//		fileName has been set
//  Postconditions:
//		sequence - populated with sequence from file
void MultipleAlignmentFile::populate() {

	ifstream inputFile( fileName);
	string line;

	getline(inputFile, line);
	vector<string> firstLineTokens;
	StringUtilities::split(line, ':', firstLineTokens);
	vector<string> positions;
	StringUtilities::split(firstLineTokens.back(), '-', positions);
	startPosition = atoi(positions.front().c_str());
		
	while(getline(inputFile, line)) {
		if (line.empty())
			continue;
		vector<string> tokens;
		StringUtilities::split(line, '\t', tokens);	
		if ("hg18" == tokens.front()) {
			string humanSeq = tokens.back();

			// Get the dog sequence
			getline(inputFile, line);
			vector<string> dogTokens;
			StringUtilities::split(line, '\t', dogTokens);	
			string  dogSeq = dogTokens.back();

			// Get the mouse sequence
			getline(inputFile, line);
			vector<string> mouseTokens;
			StringUtilities::split(line, '\t', mouseTokens);	
			string  mouseSeq = mouseTokens.back();

			// Create the residue for the sequences
			for (size_t i=0; i < humanSeq.length(); i++) {
				stringstream ss;
				ss
					<< humanSeq.at(i)
					<< dogSeq.at(i)
					<< mouseSeq.at(i);
				sequence.push_back(ss.str());
			}
		}
	}

	inputFile.close();
}
