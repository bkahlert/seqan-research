// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013


void write_to_file(Match_found & seed_found, StringSet<String<char>> & proteinID, StringSet<String<char>> & readID){

	ofstream outfile;
    outfile.open("output_example.txt");
    if (outfile.is_open()){
		outfile << "read_id\tbegin_read\tend_read\tprotein_id\tbegin_protein\tend_protein"<<endl;  
		for (int position=0;position<length(seed_found.position_read);++position){
			outfile << readID[seed_found.position_read[position]]<<"\t"<<seed_found.begin_read[position]<<"\t"<<
				seed_found.end_read[position]<<"\t"<<proteinID[seed_found.position_protein[position]]<<"\t"<<
				seed_found.begin_protein[position]<<"\t"<<seed_found.end_protein[position]<<endl;
		}
		outfile.close();
	}
}