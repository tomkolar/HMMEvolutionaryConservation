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

#ifndef MULTIPLEALIGNMENTFILE_H
#define MULTIPLEALIGNMENTFILE_H

#include <string>
#include <vector>
using namespace std;

class MultipleAlignmentFile {

public:

	// Constuctors
	// ==============================================
	MultipleAlignmentFile(); 
	MultipleAlignmentFile(string fileName);  

	// Destructor
	// =============================================
	virtual ~MultipleAlignmentFile();

	// Public Methods
	// =============================================

	// Public Accessors
	// =============================================
	const int getSequenceLength();  // length of dnaSequence
	const int getStartPosition();  // start position on chromosome
	string& MultipleAlignmentFile::getFileName();
	vector<string>& MultipleAlignmentFile::getSequence();

private:
	// Attributes
	// =============================================
    string fileName;
	int startPosition;
    vector<string> sequence;

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
    void populate();

};

#endif // MULTIPLEALIGNMENTFILE_H 
