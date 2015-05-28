#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

int main(int argc, char const ** argv)
{
	seqan::ArgumentParser parser("Read info");
	addArgument(parser, seqan::ArgParseArgument(
				seqan::ArgParseArgument::STRING, "TEXT"));

    std::cout << "Counting sequences!" << std::endl;

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    seqan::CharString fileName;
    getArgumentValue(fileName, parser, 0);

    std::cout << fileName << std::endl;

    return 1;
}
