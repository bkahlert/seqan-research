//lagan_functions.h

#ifndef LAGAN_FUNCTIONS_H
#define LAGAN_FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/score.h>
#include <seqan/seeds.h>
#include <seqan/align.h>

using namespace std;
using namespace seqan;

bool get_global_seed_chain(seqan::String<seqan::Seed<seqan::Simple> > &seedChain, seqan::DnaString &seq1, seqan::DnaString &seq2, std::pair<unsigned, unsigned> startPos, std::pair<unsigned, unsigned> endPos, unsigned q) ;

bool readFASTA(char const* path, seqan::CharString  &id, seqan::DnaString &seq);

std::string get_file_contents(const char* filepath) ;

bool writeSeedPositions(std::vector< std::pair<unsigned, unsigned> > &array1, std::vector< std::pair<unsigned, unsigned> > &array2, const char* filepath) ;

#endif
