// ==========================================================================
//                                prototyp_1_0
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Your Name <your.email@example.net>
// ==========================================================================

#include <iostream>
#include <fstream>
#include <algorithm>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
using namespace std;
using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Typedef	======	TAbsoluteMatrix	TFrequencyMatrix
// --------------------------------------------------------------------------

typedef String<ProfileChar<Dna> > TAbsoluteMatrix;
typedef String<String<double> > TFrequencyMatrix;

// --------------------------------------------------------------------------
// Class JasparRecord
// --------------------------------------------------------------------------

struct JasparRecord {
	CharString id;
	CharString name;
	TAbsoluteMatrix matrix;
	JasparRecord() {
	}
};

// --------------------------------------------------------------------------
// Class Variant
// --------------------------------------------------------------------------

struct Variant {
	unsigned position;
	Dna newNucleotide;
	Variant() {
	}
};

// --------------------------------------------------------------------------
// Class FoundBindingSite
// --------------------------------------------------------------------------

struct FoundBindingSite {
	unsigned long position;
	double referenceScore;
	double patientScore;
	FoundBindingSite() {
	}
};

// --------------------------------------------------------------------------
// Tags	======	Jaspar	ProductScore	MatchScore	AddPseudoCount
// --------------------------------------------------------------------------

struct Jaspar_;
typedef Tag<Jaspar_> Jaspar;

struct ProductScore_;
typedef Tag<ProductScore_> ProductScore;

struct MatchScore_;
typedef Tag<MatchScore_> MatchScore;

struct AddPseudoCount_;
typedef Tag<AddPseudoCount_> AddPseudoCount;

// --------------------------------------------------------------------------
// Typedef	======	TSequence	TVariantList	TFoundBindingSiteList	TSequenceList
// --------------------------------------------------------------------------

typedef String<Dna> TSequence;
typedef String<Variant> TVariantList;
typedef String<FoundBindingSite> TFoundBindingSiteList;
typedef String<TSequence> TSequenceList;

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function clear()
// Overload for JapsarRecord
// --------------------------------------------------------------------------

void clear(JasparRecord & record) {
	clear(record.id);
	clear(record.name);
	clear(record.matrix);
}

// --------------------------------------------------------------------------
// Function readRecord()
// Overload for JapsarRecord
// --------------------------------------------------------------------------

template<typename TStream, typename TPass>
int readRecord(JasparRecord & record, RecordReader<TStream, TPass> & reader,
		Jaspar const & /*tag*/) {
	clear(record);
	int res = 0;
	// Jaspar records look like this:
	//
	// ><id> <name>
	// A  [ x x x ... x ]
	// C  [ x x x ... xx ]
	// G  [xx x xx ... x ]
	// T  [xx xx xx ... x ]

	// skip until record begin
	res = skipUntilChar(reader, '>');
	if (res != 0)
		return res;
	// <id>
	if (goNext(reader)) // skip the '>'
		return EOF_BEFORE_SUCCESS;
	res = readUntilChar(record.id, reader, ' ');
	if (res != 0)
		return res;
	if (goNext(reader))
		return EOF_BEFORE_SUCCESS;
	// <name>
	res = readUntilChar(record.name, reader, '\n');
	if (res != 0)
		return res;
	if (goNext(reader))
		return EOF_BEFORE_SUCCESS;

	for (int rowIndex = 0; rowIndex < 4; ++rowIndex) { // for 0 to 3 meaning for A, C, G, T
		res = skipUntilChar(reader, '[');
		if (res != 0)
			return res;
		if (goNext(reader)) // skip the '['
			return EOF_BEFORE_SUCCESS;
		res = skipWhitespaces(reader);
		if (res != 0)
			return res;
		int columnIndex = 0;
		while (value(reader) != ']') {
			CharString bufferString;
			res = readUntilWhitespace(bufferString, reader);
			if (res != 0)
				return res;
			if (rowIndex == 0)
				resize(record.matrix, length(record.matrix) + 1);
			if (!lexicalCast2<unsigned> (
					record.matrix[columnIndex].count[rowIndex], bufferString))
				return 1; // Could not cast or could not write into matrix!
			res = skipWhitespaces(reader);
			if (res != 0)
				return res;
			++columnIndex;
		}
	}

	return 0;
}

// --------------------------------------------------------------------------
// Function absoluteToFrequencyMatrix()
// --------------------------------------------------------------------------

TFrequencyMatrix absoluteToFrequencyMatrix(TAbsoluteMatrix matrix) {
	TFrequencyMatrix myFrequencyMatrix;
	resize(myFrequencyMatrix,length(matrix));
	unsigned columnSum = matrix[0].count[0] + matrix[0].count[1] + matrix[0].count[2] + matrix[0].count[3];
	for (unsigned columnIndex=0; columnIndex<length(matrix);++columnIndex) {
		resize(myFrequencyMatrix[columnIndex],4);
		myFrequencyMatrix[columnIndex][0] = (double) matrix[columnIndex].count[0] / columnSum;
		myFrequencyMatrix[columnIndex][1] = (double) matrix[columnIndex].count[1] / columnSum;
		myFrequencyMatrix[columnIndex][2] = (double) matrix[columnIndex].count[2] / columnSum;
		myFrequencyMatrix[columnIndex][3] = static_cast<double>( matrix[columnIndex].count[3] ) / columnSum;
	}
	return myFrequencyMatrix;
}

TFrequencyMatrix absoluteToFrequencyMatrix(TAbsoluteMatrix matrix, AddPseudoCount const & /* tag */) {
	TFrequencyMatrix myFrequencyMatrix;
	resize(myFrequencyMatrix,length(matrix));
	unsigned columnSum = matrix[0].count[0] + matrix[0].count[1] + matrix[0].count[2] + matrix[0].count[3] + 4;
	for (unsigned columnIndex=0; columnIndex<length(matrix);++columnIndex) {
		resize(myFrequencyMatrix[columnIndex],4);
		myFrequencyMatrix[columnIndex][0] = (double) (matrix[columnIndex].count[0]+1) / columnSum;
		myFrequencyMatrix[columnIndex][1] = (double) (matrix[columnIndex].count[1]+1) / columnSum;
		myFrequencyMatrix[columnIndex][2] = (double) (matrix[columnIndex].count[2]+1) / columnSum;
		myFrequencyMatrix[columnIndex][3] = (double) (matrix[columnIndex].count[3]+1) / columnSum;
	}
	return myFrequencyMatrix;
}
// --------------------------------------------------------------------------
// Function createFrequencyMatrix()
// --------------------------------------------------------------------------

TAbsoluteMatrix createAbsoluteMatrix() {
	fstream stream("jasparMatrix_prototyp_1_0.txt", ios::binary | ios::in);
	// Read file.
	RecordReader < std::fstream, SinglePass<> > reader(stream);
	JasparRecord myRecord;
	readRecord(myRecord, reader, Jaspar());
	return myRecord.matrix;
}

// --------------------------------------------------------------------------
// Function createSequence()
// --------------------------------------------------------------------------

TSequence createSequence() {
	TSequence mySequence;
	appendValue(mySequence, 'A'); // 0
	appendValue(mySequence, 'G');
	appendValue(mySequence, 'T');
	appendValue(mySequence, 'C');
	appendValue(mySequence, 'C');
	appendValue(mySequence, 'T');
	appendValue(mySequence, 'A');
	appendValue(mySequence, 'A');
	appendValue(mySequence, 'T');
	appendValue(mySequence, 'T');
	appendValue(mySequence, 'T'); // 10
	appendValue(mySequence, 'G');
	appendValue(mySequence, 'G');
	appendValue(mySequence, 'C');
	appendValue(mySequence, 'C');
	appendValue(mySequence, 'C'); // 15
	appendValue(mySequence, 'A');
	appendValue(mySequence, 'A');
	appendValue(mySequence, 'A');
	appendValue(mySequence, 'T');
	appendValue(mySequence, 'T'); // 20
	appendValue(mySequence, 'T');
	appendValue(mySequence, 'G');
	appendValue(mySequence, 'G');
	appendValue(mySequence, 'G');
	appendValue(mySequence, 'G'); // 25
	appendValue(mySequence, 'G');
	appendValue(mySequence, 'G');
	return mySequence;
}

// --------------------------------------------------------------------------
// Function createVariantList()
// --------------------------------------------------------------------------

TVariantList createVariantList() {
	Variant myVariant1;
	myVariant1.position = 5;
	myVariant1.newNucleotide = 'A';
	Variant myVariant2;
	myVariant2.position = 13;
	myVariant2.newNucleotide = 'A';
	Variant myVariant3;
	myVariant3.position = 27;
	myVariant3.newNucleotide = 'C';
	TVariantList myVariantList;
	appendValue(myVariantList, myVariant1);
	appendValue(myVariantList, myVariant2);
	appendValue(myVariantList, myVariant3);
	return myVariantList;
}

// --------------------------------------------------------------------------
// Function getProductScore()
// --------------------------------------------------------------------------

double getProductScore(TSequence sequence, TFrequencyMatrix matrix,
		unsigned sequencePosition) {
	double myCurrentScore = 1;
	for (unsigned matrixColumn = 0; matrixColumn < length(matrix); ++matrixColumn) {
		if (sequence[sequencePosition + matrixColumn] == 'A')
			myCurrentScore = myCurrentScore * matrix[matrixColumn][0];
		else if (sequence[sequencePosition + matrixColumn] == 'C')
			myCurrentScore = myCurrentScore * matrix[matrixColumn][1];
		else if (sequence[sequencePosition + matrixColumn] == 'G')
			myCurrentScore = myCurrentScore * matrix[matrixColumn][2];
		else if (sequence[sequencePosition + matrixColumn] == 'T')
			myCurrentScore = myCurrentScore * matrix[matrixColumn][3];
	}

	return myCurrentScore;
}

// --------------------------------------------------------------------------
// Function getMatchScore()
// --------------------------------------------------------------------------

double getMatchScore(TSequence sequence, TAbsoluteMatrix matrix,
		unsigned sequencePosition) {
	double myCurrentScore = 0;
	double currentMin = 0;
	double currentMax = 0;
	for (unsigned matrixColumn = 0; matrixColumn < length(matrix); ++matrixColumn) {
		if (sequence[sequencePosition + matrixColumn] == 'A')
			myCurrentScore += matrix[matrixColumn].count[0];
		else if (sequence[sequencePosition + matrixColumn] == 'C')
			myCurrentScore += matrix[matrixColumn].count[1];
		else if (sequence[sequencePosition + matrixColumn] == 'G')
			myCurrentScore += matrix[matrixColumn].count[2];
		else if (sequence[sequencePosition + matrixColumn] == 'T')
			myCurrentScore += matrix[matrixColumn].count[3];

		currentMin += *min(begin(matrix[matrixColumn].count),end(matrix[matrixColumn].count));
		currentMax += *max(begin(matrix[matrixColumn].count),end(matrix[matrixColumn].count));
	}
	return ((myCurrentScore - currentMin) / (currentMax - currentMin));
}

// --------------------------------------------------------------------------
// Function getAllPossibleSequencesOfGivenLength()
// --------------------------------------------------------------------------

// Helper function for getAllPossibleSequencesOfGivenLength().
void getAllPossibleSequencesOfGivenLengthHelper(TSequenceList & mySequenceList,
		TSequence &current, unsigned pos) {
	if (pos < length(current)) {
		current[pos] = 'A';
		getAllPossibleSequencesOfGivenLengthHelper(mySequenceList, current, pos + 1);
		current[pos] = 'C';
		getAllPossibleSequencesOfGivenLengthHelper(mySequenceList, current, pos + 1);
		current[pos] = 'G';
		getAllPossibleSequencesOfGivenLengthHelper(mySequenceList, current, pos + 1);
		current[pos] = 'T';
		getAllPossibleSequencesOfGivenLengthHelper(mySequenceList, current, pos + 1);
	} else {
		appendValue(mySequenceList, current);
	}
}

TSequenceList getAllPossibleSequencesOfGivenLength(unsigned length) {
	TSequence current;
	TSequenceList result;
	resize(current, length);
	getAllPossibleSequencesOfGivenLengthHelper(result, current, 0);
	return result;
}

// --------------------------------------------------------------------------
// Function determineThreshold()
// --------------------------------------------------------------------------

double determineThreshold(TFrequencyMatrix matrix) {
	TSequenceList allPossibleSequences;
	allPossibleSequences = getAllPossibleSequencesOfGivenLength(length(matrix));
	String<double> possibleScores;

	typedef Iterator<TSequenceList>::Type TSequenceListIterator;
	for (TSequenceListIterator it = begin(allPossibleSequences); it != end(
			allPossibleSequences); goNext(it)) {
		appendValue(possibleScores, getProductScore(value(it), matrix, 0));
	}

	typedef Iterator<double>::Type TDoubleIterator;
	TDoubleIterator beginIt = begin(possibleScores);
	TDoubleIterator endIt = end(possibleScores);
	sort(beginIt, endIt);
	beginIt = begin(possibleScores);
	endIt = end(possibleScores);
	double scoreSum=0;
	while (endIt != beginIt) {
		goPrevious(endIt);
		scoreSum += value(endIt);
		if (scoreSum >= 0.95)
			break;
	}
	return value(endIt);
}

// --------------------------------------------------------------------------
// Function searchRelevantBindingsites()
// --------------------------------------------------------------------------

TFoundBindingSiteList searchRelevantBindingsites(TFrequencyMatrix matrix,
		TSequence referenceSequence, TVariantList variantList,
		double threshold, ProductScore const & /* tag*/) {
	TFoundBindingSiteList myFoundBindingSiteList;
	FoundBindingSite myCurrentBindingSite;
	TSequence patientSequence = referenceSequence;
	for (unsigned variantListPosition = 0; variantListPosition < length(
			variantList); ++variantListPosition) {
		patientSequence[variantList[variantListPosition].position]
				= variantList[variantListPosition].newNucleotide;
	}
	for (unsigned sequencePosition = 0; sequencePosition < length(referenceSequence)
			- length(matrix) + 1; ++sequencePosition) {
		double currentReferenceScore = getProductScore(referenceSequence,
				matrix, sequencePosition);
		double currentPatientScore = getProductScore(patientSequence, matrix,
				sequencePosition);
		if (((currentReferenceScore >= threshold) or (currentPatientScore >= threshold))
				and (currentReferenceScore != currentPatientScore)) {
			myCurrentBindingSite.position = sequencePosition;
			myCurrentBindingSite.referenceScore = currentReferenceScore;
			myCurrentBindingSite.patientScore = currentPatientScore;
			appendValue(myFoundBindingSiteList, myCurrentBindingSite);
		}
	}
	return myFoundBindingSiteList;
}

TFoundBindingSiteList searchRelevantBindingsites(TAbsoluteMatrix matrix,
		TSequence referenceSequence, TVariantList variantList, double threshold,
		MatchScore const & /* tag*/) {
	TFoundBindingSiteList myFoundBindingSiteList;
	FoundBindingSite myCurrentBindingSite;
	TSequence patientSequence = referenceSequence;
	for (unsigned variantListPosition = 0; variantListPosition < length(
			variantList); ++variantListPosition) {
		patientSequence[variantList[variantListPosition].position]
				= variantList[variantListPosition].newNucleotide;
	}
	for (unsigned sequencePosition = 0; sequencePosition < length(referenceSequence)
			- length(matrix) + 1; ++sequencePosition) {
		double currentReferenceScore = getMatchScore(referenceSequence,
				matrix, sequencePosition);
		double currentPatientScore = getMatchScore(patientSequence, matrix,
				sequencePosition);
		if (((currentReferenceScore >= threshold) or (currentPatientScore >= threshold))
				and (currentReferenceScore != currentPatientScore)) {
			myCurrentBindingSite.position = sequencePosition;
			myCurrentBindingSite.referenceScore = currentReferenceScore;
			myCurrentBindingSite.patientScore = currentPatientScore;
			appendValue(myFoundBindingSiteList, myCurrentBindingSite);
		}
	}
	return myFoundBindingSiteList;
}

// --------------------------------------------------------------------------
// Function printFoundBindingSiteList()
// --------------------------------------------------------------------------

void printFoundBindingSiteList(TFoundBindingSiteList myFoundBindingSiteList, double myThreshold) {
	cout << "=== TransVar =====================================" << endl;
	cout << endl;
	cout << "Threshold = " << myThreshold << endl;
	cout << endl;
	cout << "pos" << "\t" << "refScore" << "\t" << "patScore" << "\t"
			<< "dif_abs " << "\t" << "dif_rel " << endl;
	for (unsigned foundBindingSiteIndex = 0; foundBindingSiteIndex < length(
			myFoundBindingSiteList); ++foundBindingSiteIndex) {
		cout << myFoundBindingSiteList[foundBindingSiteIndex].position << "\t"
				<< myFoundBindingSiteList[foundBindingSiteIndex].referenceScore
				<< "\t"
				<< myFoundBindingSiteList[foundBindingSiteIndex].patientScore
				<< "\t"
				<< myFoundBindingSiteList[foundBindingSiteIndex].patientScore
						- myFoundBindingSiteList[foundBindingSiteIndex].referenceScore
				<< "\t"
				<< (myFoundBindingSiteList[foundBindingSiteIndex].patientScore
						- myFoundBindingSiteList[foundBindingSiteIndex].referenceScore)
						/ myFoundBindingSiteList[foundBindingSiteIndex].referenceScore
				<< endl;
	}
	cout << endl;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main() {
	TAbsoluteMatrix myAbsoluteMatrix = createAbsoluteMatrix();
	TFrequencyMatrix myFrequencyMatrix = absoluteToFrequencyMatrix(myAbsoluteMatrix,AddPseudoCount());
	TSequence mySequence = createSequence();
	TVariantList myVariantList = createVariantList();
//	double myThreshold = determineThreshold(myFrequencyMatrix);
	double myThreshold = 7.34467e-07;
	printFoundBindingSiteList(searchRelevantBindingsites(myFrequencyMatrix, mySequence, myVariantList, myThreshold, ProductScore()),myThreshold);
	printFoundBindingSiteList(searchRelevantBindingsites(myAbsoluteMatrix, mySequence, myVariantList, 0.5, MatchScore()),0.5);
	return 0;
}

