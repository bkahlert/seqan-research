
#ifndef PEMER_LITE_H_
#define PEMER_LITE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <algorithm>
#include <iterator>
#include <iostream>
#include <seqan/bam_io.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <math.h>
#include <vector>
#include <fstream>

struct candidate;
struct modifier;
int median(candidate &input);
unsigned average(candidate &input);
seqan::ArgumentParser::ParseResult commands(modifier &options, int argc, char const ** argv);
int analyse(candidate &input);
int saminput(candidate &save, char *file);
int devide(candidate &input,candidate &result,modifier &options);
int tsv(std::vector<std::vector<unsigned>> &input,candidate &ref, char *out);
int findCluster(candidate &input, std::vector<std::vector<unsigned> > &dest, candidate::form indel);

#endif  // PEMER_LITE_HEADER_H_
