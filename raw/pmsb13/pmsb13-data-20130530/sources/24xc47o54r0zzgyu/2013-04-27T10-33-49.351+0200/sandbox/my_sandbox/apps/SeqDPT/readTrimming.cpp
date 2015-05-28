#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/seq_io.h>
#include <seqan/file.h>      // to stream a CharString into cout
#include "readTrimming.h"

// Trimming methods
// TODO: add support for minimum length and removing smaller sequences (or replace with one N).
// ----------------------------------------------------------------------------
// Simple trimming mechanism: Cut off bases as long as they are below our quality.
template <typename TSeq, typename TQual>
int _trimRead(TSeq& seq, TQual& qual, int const cutoff, int const min_length, Standard const &)
{
	int cut_pos = length(seq);
	int ret;
	for (int i=length(seq)-1; i >= min_length; i--)
		if ((seqan::ordValue(qual[i]) - 33) < cutoff)
			cut_pos = i;

	return cut_pos;
}

template <typename TSeq>
int _trimRead(TSeq& seq, int const cutoff, int const min_length, Standard const &)
{
	int cut_pos = length(seq);
	int ret;
	for (int i=length(seq)-1; i >= min_length; i--)
		if (seqan::getQualityValue(seq[i]) < cutoff)
			cut_pos = i;

	return cut_pos;
}

// Trimming mechanism used in BWA. Trim to argmax_x sum_{i=x+1}^l {cutoff - q_i}
template <typename TSeq, typename TQual>
int _trimRead(TSeq& seq, TQual& qual, int const cutoff, int const min_length, BWA const &)
{
	int max_arg = length(seq), sum = 0, max = 0;
	for (int i=length(seq)-1; i >= min_length; i--)
	{
		sum += cutoff - (seqan::ordValue(qual[i]) - 33);
		if (sum < 0)
			break;
		if (sum > max)
		{
			max = sum;
			max_arg = i;
		}
	}

	return max_arg;
}

template <typename TSeq>
int _trimRead(TSeq& seq, int const cutoff, int const min_length, BWA const &)
{
	int max_arg = length(seq), sum = 0, max = 0;
	for (int i=length(seq)-1; i >= min_length; i--)
	{
		sum += cutoff - seqan::getQualityValue(seq[i]);
		if (sum < 0)
			break;
		if (sum > max)
		{
			max = sum;
			max_arg = i;
		}
	}

	return max_arg;
}

// Trim by shifting a window over the text and calculating mean/median in it.
template <typename TSeq, typename TQual>
int _trimRead(TSeq& seq, TQual& qual, int const cutoff, int const min_length, Median const & spec)
{
	int window = spec.window;
	int bad_pos = length(seq);

	// Shift window over read and keep mean quality.
	for (int i=min_length; i < length(seq)-window; i++)
	{
		double avg = 0;
		bad_pos = length(seq);
		// Calculate mean quality in the window.
		for (int j=0; j < window; j++)
		{
			avg += seqan::ordValue(qual[i+j]);
			// Save position of first bad base in the window.
			// (We start to cut here, if it turns out the whole window is bad.)
			if (bad_pos == length(seq) && seqan::ordValue(qual[i+j] - 33) < cutoff)
				bad_pos = i+j;
		}
		avg = (avg/window) - 33;

		// If this window was bad, we stop here and cut the tail.
		if (avg < cutoff)
			break;
	}

	// We will cut if we found a bad window <=> if bad_pos was changed.
	return bad_pos;
}

template <typename TSeq>
int _trimRead(TSeq& seq, int const cutoff, int const min_length, Median const & spec)
{
	/*int window = spec.window;
	int bad_pos = length(seq);

	// Shift window over read and keep mean quality.
	for (int i=min_length; i < length(seq)-window; i++)
	{
		double avg = 0;
		bad_pos = length(seq);
		// Calculate mean quality in the window.
		for (int j=0; j < window; j++)
		{
			avg += seqan::getQualityValue(seq[i+j]);
			// Save position of first bad base in the window.
			// (We start to cut here, if it turns out the whole window is bad.)
			if (bad_pos == length(seq) && seqan::getQualityValue(seq[i+j]) < cutoff)
				bad_pos = i+j;
		}
		avg = (avg/window);

		// If this window was bad, we stop here and cut the tail.
		if (avg < cutoff)
			break;
	}

	// We will cut if we found a bad window <=> if bad_pos was changed.
	return bad_pos;*/

	int window = spec.window;
	int bad_pos = min_length;
	double avg = 0;

	for (int i=0; i < window; i++)
		avg += seqan::getQualityValue(seq[i]);
	avg /= window;

	for (int i=window; i < length(seq) && avg >= cutoff; i++)
	{
		int drop = i-window;
		avg = avg + (seqan::getQualityValue(seq[i]) - seqan::getQualityValue(seq[drop]))/window;
		if (bad_pos == drop)
		{
			bad_pos++;
			while (bad_pos < length(seq) && seqan::getQualityValue(seq[bad_pos]) >= cutoff)
				bad_pos++;
		}
	}

	return bad_pos;
}

// Trimming functions. With or without dedicated quality string.
// Return the number of bases trimmed off the end.
template <typename TSeq, typename TQual, typename TSpec>
int trimRead(TSeq& seq, TQual& qual, int const cutoff, int const min_length, TSpec const & spec)
{
	int ret, cut_pos;
	cut_pos = _trimRead(seq, qual, cutoff, min_length, spec);

	ret = length(seq) - cut_pos;
	seqan::erase(seq , cut_pos, length(seq));
	seqan::erase(qual, cut_pos, length(qual));

	return ret;
}

template <typename TSeq, typename TSpec>
int trimRead(TSeq& seq, int const cutoff, int const min_length, TSpec const & spec)
{
	int ret, cut_pos;
	cut_pos = _trimRead(seq, cutoff, min_length, spec);

	ret = length(seq) - cut_pos;
	seqan::erase(seq, cut_pos, length(seq));

	return ret;
}

// trimBatch applies the trimRead function to each Sequence in a StringSet
// returns the number of trimmed sequences.
template <typename TSeq, typename TQual, typename TSpec>
int trimBatch(seqan::StringSet<TSeq> seqSet, seqan::StringSet<TQual> qualSet, int const cutoff, int const min_length, TSpec const & spec)
{
	int trimmedReads = 0;
	for (int i=0; i < length(seqSet); i++)
	{
		if (trimRead(value(seqSet, i), value(qualSet, i), cutoff, min_length, spec) > 0)
			trimmedReads++;
	}
	return trimmedReads;
}

template <typename TSeq, typename TSpec>
int trimBatch(seqan::StringSet<TSeq> seqSet, int const cutoff, int const min_length, TSpec const & spec)
{
	int trimmedReads = 0;
	for (int i=0; i < length(seqSet); i++)
	{
		if (trimRead(value(seqSet, i), cutoff, min_length, spec) > 0)
			trimmedReads++;
	}
	return trimmedReads;
}


// Trims each read in a fastq file. Will probably be removed.
template <typename TSpec>
int trimFile(seqan::String<char> const inFile, seqan::String<char> const outFile,
				int const cutoff, int const min_length, TSpec const & spec)
{
	seqan::SequenceStream inStream(toCString(inFile));
    if (!isGood(inStream))
    {
        std::cerr << "ERROR: Could not open the file" << inFile << "\n";
        return 1;
    }

    //seqan::SequenceStream outStream(toCString(outFile), seqan::SequenceStream::WRITE);
   /* seqan::SequenceStream outStream("trimmed.fq", seqan::SequenceStream::WRITE);

    if (!isGood(outStream))
    {
		std::cerr << "ERROR: Could not open the file" << outFile << "\n";
		return 1;
    }*/

    seqan::StringSet<seqan::String<char> > ids;
    seqan::StringSet<seqan::String<seqan::Dna5Q> > seqs;
    typedef seqan::Iterator<seqan::StringSet<seqan::String<seqan::Dna5Q> > > TDna5QSet;
    int count = 0;
    while (!atEnd(inStream))
    {
    	if (readBatch(ids, seqs, inStream, 1000) == 0)
    	{
    		for (int i=0; i < length(seqs); i++)
    		{
    			//seqan::String<seqan::Dna5Q> tmp;
    			//resize(tmp, length(value(seqs,i)));
    			//assign(tmp, value(seqs,i));
    			if (trimRead(value(seqs, i), cutoff, min_length, spec) > 0){
    				//std::cout << value(ids,i) << "\n" << tmp << "\n" << value(seqs,i) << "\n";
    				count++;
    			}
    		}

    		/*if (writeAll(outStream, ids, seqs) != 0)
    		{
    	        std::cerr << "ERROR: Could not write to file!\n";
    	        return 1;
    		}*/
    	}
    }

    std::cout << "Getrimmte Reads: " << count << std::endl;
}
