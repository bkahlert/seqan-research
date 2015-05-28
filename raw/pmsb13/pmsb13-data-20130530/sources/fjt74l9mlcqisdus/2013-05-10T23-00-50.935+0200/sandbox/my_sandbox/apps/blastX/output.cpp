// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013
#include "own_functions.h"

void write_to_file(Match_found & seed_found, StringSet<String<char> > & proteinID, StringSet<String<char> > & readID){
	
	ofstream outfile;
	outfile.open("output_example.txt");
    if (outfile.is_open()){
		outfile << "read_id\tbegin_read\tend_read\tprotein_id\tbegin_protein\tend_protein\tframe\tside"<<endl;  
		for (int position=0;position<length(seed_found.position_read);++position){
			int read = seed_found.position_read[position]/6;
			int protein = seed_found.position_protein[position];
			int frame = seed_found.position_read[position]%6;
			int side = frame/3;
			outfile << readID[read]
			<<"\t"<<seed_found.begin_read[position]
			<<"\t"<<seed_found.end_read[position]
			<<"\t"<<proteinID[protein]
			<<"\t"<<seed_found.begin_protein[position]
			<<"\t"<<seed_found.end_protein[position]
			<<"\t"<<frame
			<<"\t"<<side<<endl;

		}
		outfile.close();
	}
	
}


void write_to_file(Match_found & verify_found, StringSet<String<char> > & proteinID, StringSet<String<char> > & readID,StringSet<String<Dna>> & Reads){
	
	ofstream outfile;
	outfile.open("output_example_verify.txt");
    if (outfile.is_open()){
		outfile << "read_id\tread_sequence\tprotein_id\tbegin_protein\tend_protein\tframe\tside"<<endl;  
		for (int position=0;position<length(verify_found.position_read);++position){
			int read = verify_found.position_read[position]/6;
			int protein = verify_found.position_protein[position];
			int frame = verify_found.position_read[position]%6;
			int side = frame/3;
			outfile << readID[read]
			<<"\t"<<Reads[read]
			<<"\t"<<proteinID[protein]
			<<"\t"<<verify.begin_protein[position]
			<<"\t"<<verify.end_protein[position]
			<<"\t"<<frame
			<<"\t"<<side<<endl;

		}
		outfile.close();
	}
	
}
