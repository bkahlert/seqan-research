// BLASTX IMPLEMENTIERUNG VON ANNKATRIN BRESSIN UND MARJAN FAIZI
// SOFTWAREPROJEKT VOM 2.4. - 29.5.2012
// VERSION VOM 04.MAI.2013

/**@brief output.cpp
*Beinhaltet Funktionen, die die Ergebnisse in eine Textdatei schreiben
*/

#include "own_functions.h"


void write_to_file(Match_found & seed_found, StringSet<String<char> > & proteinID, StringSet<String<char> > & readID,String<unsigned int> numb_alp,int alp){
	/**@brief write_to_file schreibt die gefundenen matches in eine Textdatei
	*@param seed_found Instanz der Klasse Match_found enthaelt die Informationen der gefundenen matches 
	*@param proteinID IDs der einzelnen Proteinsequenzen
	*@param readID IDs der einzelnen Reads
	*@param numb_alp Speichert die Zahlen der ausgesuchten Alphabete
	*@param alp Aktuelle Position des benutzten Alphabets
	*@return Hat kein Rueckgabewert
	*/
	/*
	StringSet<String<char> > matrix_name;
	appendValue(matrix_name,"Blosum30");
	appendValue(matrix_name,"Blosum45");
	appendValue(matrix_name,"Blosum80");
	appendValue(matrix_name,"Pam120");
	appendValue(matrix_name,"Pam200");
	appendValue(matrix_name,"Pam250");
	appendValue(matrix_name,"Pam40");
	appendValue(matrix_name,"Vtml200");
	*/
	string matrix[8] = {"Blosum30","Blosum45","Blosum80","Pam120","Pam200","Pam250","Pam40","Vtml200"};
	
	string file = "output" + matrix[numb_alp[alp]-49] + ".txt";
	ofstream outfile;
	outfile.open(toCString(file));
    if (outfile.is_open()){
		outfile << "read_id\tbegin_read\tend_read\tprotein_id\tbegin_protein\tend_protein\tframe\tside"<<endl;  
		for (int position=0;position<length(seed_found.position_read);++position){
			int read = (seed_found.position_read[position]/6);
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


void write_to_file(Match_found & verify_found, StringSet<String<char> > & proteinID, StringSet<String<char> > & readID,StringSet<String<Dna> > & reads,String<unsigned int> numb_alp,int alp){
	/**@brief write_to_file schreibt die verifizierten matches in eine Textdatei
	*@param verify_found Instanz der Klasse Match_found enthaelt die Informationen der verifizierten matches 
	*@param proteinID IDs der einzelnen Proteinsequenzen
	*@param readID IDs der einzelnen Reads
	*@param reads Original Reads(nicht uebersetzt)
	*@param numb_alp Speichert die Zahlen der ausgesuchten Alphabete
	*@param alp Aktuelle Position des benutzten Alphabets
	*@return Hat kein Rueckgabewert
	*/
	string matrix[8] = {"Blosum30","Blosum45","Blosum80","Pam120","Pam200","Pam250","Pam40","Vtml200"};
	string file = "output" + matrix[numb_alp[alp]-49] + "verify.txt";
	ofstream outfile;
	outfile.open(toCString(file));
    if (outfile.is_open()){
		outfile << "read_id\tread_sequence\tprotein_id\tbegin_protein\tend_protein\tframe\tside"<<endl;  
		for (int position=0;position<length(verify_found.position_read);++position){
			int read = verify_found.position_read[position]/6;
			int protein = verify_found.position_protein[position];
			int frame = verify_found.position_read[position]%6;
			int side = frame/3;
			outfile << readID[read]
			<<"\t"<<reads[read]
			<<"\t"<<proteinID[protein]
			<<"\t"<<verify_found.begin_protein[position]
			<<"\t"<<verify_found.end_protein[position]
			<<"\t"<<frame
			<<"\t"<<side<<endl;
		}
		outfile.close();
	}
}
