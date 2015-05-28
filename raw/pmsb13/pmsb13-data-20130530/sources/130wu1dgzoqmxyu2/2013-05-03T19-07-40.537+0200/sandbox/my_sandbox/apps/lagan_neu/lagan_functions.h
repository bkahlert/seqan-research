#ifndef LAGAN_FUNCTIONS_H
#define LAGAN_FUNCTIONS_H

bool get_global_seed_chain(seqan::String<seqan::Seed<seqan::Simple> > &seedChain, seqan::DnaString &seq1, seqan::DnaString &seq2, std::pair<unsigned, unsigned> startPos, std::pair<unsigned, unsigned> endPos, unsigned q) ;

bool readFASTA(char const * path, seqan::CharString  &id, seqan::DnaString &seq) ;

std::string get_file_contents(const char* filepath) ;

bool writeSeedPositions(std::vector< std::pair<unsigned, unsigned> > &array1, std::vector< std::pair<unsigned, unsigned> > &array2, const char* filepath) ;

#endif
