/*
 *   Enums.h
 *
 *   Author: mat and jtr
 */

#ifndef FLEXBAR_ENUMS_H_
#define FLEXBAR_ENUMS_H_


/* These enums are used by almost every class. */

namespace flexbar{
	
	const unsigned int MAX_READLENGTH = 2048;
	
	enum LogLevel {
		NONE,
		ALL,
		TAB,
		MOD
	};
	
	enum TrimEnd {
		ANY,
		LEFT,
		RIGHT,
		LEFT_TAIL,
		RIGHT_TAIL
	};
	
	enum FileFormat {
		FASTA,
		FASTQ,
		CSFASTA,
		CSFASTQ
	};
	
	enum QualityType {
		SANGER,
		SOLEXA,
		ILLUMINA13
	};
	
	enum BarcodeDetect {
		BARCODE_READ,
		WITHIN_READ,
		WITHIN_READ_REMOVAL,
		BOFF
	};
	
	enum AdapterRemoval {
		NORMAL,
		AOFF
	};
	
	enum RunType {
		SINGLE,
		PAIRED,
		SINGLE_BARCODED,
		PAIRED_BARCODED
	};
	
	std::string toFormatString(FileFormat format){
		switch(format){
			case FASTA:   return ".fasta";
			case FASTQ:   return ".fastq";
			case CSFASTA: return ".csfasta";
			case CSFASTQ: return ".csfastq";
		}
		return ".unknown";
	}
}

#endif /* FLEXBAR_ENUMS_H_ */
