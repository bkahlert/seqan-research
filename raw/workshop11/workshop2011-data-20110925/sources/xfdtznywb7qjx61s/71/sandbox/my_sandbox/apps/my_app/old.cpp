/*
 *  old.cpp
 *  seqan
 *
 *  Created by Roland Krause on 13.09.11.
 *  Copyright 2011 MPI for Molecular Genetics. All rights reserved.
 *
 */

#include "old.h"

countOneMers("mississipi");

typedef String<AminoAcid> TAminoAcidString;
TAminoAcidString str = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
Iterator<String<AminoAcid> >::Type it = begin(str);
Iterator<String<AminoAcid> >::Type itEnd = end(str);

//	std::cout << str;
while(it != itEnd)
{
	if(value(it) == 'R')
	{
		std::cout<< 'A';
		} else {
			std::cout << *it;
		}
		++it;
		
		}
		
		
		showAllLetterOfMyAlphabet(AminoAcid());
		showAllLetterOfMyAlphabet(Dna());
		showAllLetterOfMyAlphabet(Dna5());
		return 0;