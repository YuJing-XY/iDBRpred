/*
 * Properties.cpp
 *
 *  Created on: Mar 2, 2016
 *      Author: nmalhis
 */

#include "Properties.h"

Properties::Properties() {
	// ----------------------- Default values
	debug = 0;
	threads = 1;
	SwissProt.assign("/data/home/nmalhis/database/DB/SP/SP");
	UniRef90.assign("/data/home/nmalhis/database/DB/uniref90filt/uniref90filt");
	inputFile.assign("input.fasta");
	outputFile.assign("output.txt");
	Local_path.assign("");
	ESpritz_path.assign("/home/nmalhis/Espritz/");
	PSIBLAST_path.assign("");
	release.assign("<Release 1.0 Dec. 15 2015>");
	// ---------------------------------------
	string line;
	ifstream fin;
	string cmd;
	string whitespaces(" \t\f\v\n\r");
	size_t pos;
	size_t last;
	set<string> labels;
	labels.insert("Input");
	labels.insert("Output");
	labels.insert("SwissProt");
	labels.insert("UniRef90");
	labels.insert("threads");
	labels.insert("ESpritz");
	labels.insert("PSIBLAST");
	labels.insert("debug");

	fin.open("MoRFchibi.properties");
	if (!fin.is_open()) {
		cerr
				<< "can\'t open MoRFchibi.properties file.\nDefault properties are used"
				<< endl;
	} else {
		getline(fin, line);
		while (!fin.eof()) {
			if (line.at(0) != '#') {
				last = line.find_last_not_of(whitespaces);
				line.erase(last + 1);
				stringstream ss(stringstream::in | stringstream::out);
				ss.str(line);
				ss >> cmd;
				if (labels.find(cmd) != labels.end()) {
					pos = line.find_last_of('\t');
					if (pos != string::npos) {
						pos++;
						while (line.at(pos) == ' ')
							pos++;
						if (pos < line.size())
							switch (cmd.at(0)) {
							case 'I':
								inputFile.assign(line.substr(pos));
								break;
							case 'O':
								outputFile.assign(line.substr(pos));
								break;
							case 'S':
								SwissProt.assign(line.substr(pos));
								break;
							case 'U':
								UniRef90.assign(line.substr(pos));
								break;
							case 'E':
								ESpritz_path.assign(line.substr(pos));
								break;
							case 'P':
								PSIBLAST_path.assign(line.substr(pos));
								break;
							case 'd':
								ss >> debug;
								break;
							case 't':
								ss >> threads;
								break;
							default:
								cout << "Unknown label";
							}
					}
				} else {
					cerr << "Error in Properties file: unknown label: <" << cmd
							<< ">" << endl;
				}
			}
			getline(fin, line);
		}
		fin.close();
		int ret = system("pwd > tmp/ld.txt");
		fin.open("tmp/ld.txt");
		getline(fin, Local_path);
		Local_path.append("/");
		fin.close();
	}
}

int Properties::clean_fasta(string &cleanOut, unsigned minSize) {
	if (minSize < 11)
		minSize = 11;
	vector<vector<char> > seqList;
	vector<string> nameList;
	vector<string> errorList;
	vector<int> errorCountList;
	int countGood = 0;
	string line;
	string cline;
	string seqError;
	vector<char> cvec;
	vector<char> seq;
	ifstream fin;
	set<char> AA;
	set<char> ignore;
	AA.insert('A');
	AA.insert('C');
	AA.insert('D');
	AA.insert('E');
	AA.insert('F');
	AA.insert('G');
	AA.insert('H');
	AA.insert('I');
	AA.insert('K');
	AA.insert('L'); // 10
	AA.insert('M');
	AA.insert('N');
	AA.insert('P');
	AA.insert('Q');
	AA.insert('R'); // 15
	AA.insert('S');
	AA.insert('T');
	AA.insert('V');
	AA.insert('W');
	AA.insert('Y');
	ignore.insert(' ');
	ignore.insert('\t');
	fin.open(inputFile.c_str());
	if (!fin.is_open()) {
	  cerr << "Fatal Error: can\'t open input file: " << inputFile << endl;
	  exit(0);
	} else {
	  getline(fin, line);
	  seq.clear();
	  while (!fin.eof()) {
	    if (line.size() > 0) {
	      if (line.at(line.size() - 1) != 13) {
		cline.assign(line);
	      } else {
		cline.assign(line.substr(0, line.size() - 1));
	      }
	      if (cline.size() > 0) {
		if (nameList.size() == 0) {
		  if (cline.at(0) != '>') {
		    cerr << "Bad fasta, line:\n\t" << cline << endl;
		    cerr << " must start with '>'" << endl;
		    exit(0);
		  }
		}
		// --------------------------------------
		if (cline.at(0) == '>') {
		  if (nameList.size() > seqList.size()) {
		    seqList.push_back(seq);
		  }
		  nameList.push_back(cline);
		  seq.clear();
		} else {
		  for (unsigned i = 0; i < cline.size(); i++) {
		    if (ignore.find(cline.at(i)) == ignore.end()) {
		      seq.push_back(cline.at(i));
		    }
		  }
		}
		// --------------------------------------
	      }
	    }
	    getline(fin, line);
	  }
	  if (nameList.size() > seqList.size()) {
	    seqList.push_back(seq);
	  }
	  fin.close();
	  int errorCount;
	  for (unsigned i = 0; i < nameList.size(); i++) {
	    errorCount = 0;
	    seqError.assign("");
	    if (seqList.at(i).size() < minSize) {
	      errorCount++;
	      stringstream ss;
	      ss << "Error 1: sequence is very short, it is just ";
	      ss << seqList.at(i).size();
	      ss << " residues, must be >= ";
	      ss << minSize;
	      ss << "\n";
	      seqError.append(ss.str());
	    }
	    for (unsigned j = 0; j < seqList.at(i).size(); j++) {
	      if (AA.find(seqList.at(i).at(j)) == AA.end()) {
		errorCount++;
		stringstream ss;
		ss << "Error ";
		ss << errorCount;
		ss << ": at index ";
		ss << (j + 1);
		ss << "; (";
		ss << seqList.at(i).at(j);
		ss << ") is not a standard amino acid\n";
		seqError.append(ss.str());
	      }
	    }
	    errorList.push_back(seqError);
	    errorCountList.push_back(errorCount);
	  }
	  errorCount = 0;
	  ofstream fout;
	  fout.open(cleanOut.c_str()); //("tmp/junk");
	  int ret = system("rm -r badSequences 2> tmp/errMsg.txt");
	  for (unsigned i = 0; i < nameList.size(); i++) {
	    if (errorCountList.at(i) == 0) {
	      fout << nameList.at(i) << endl;
	      for (unsigned j = 0; j < seqList.at(i).size(); j++) {
		fout << seqList.at(i).at(j);
	      }
	      fout << endl;
	    } else {
	      ret = system("mkdir badSequences 2> tmp/errMsg.txt");
	      errorCount++;
	      string efname;
	      ofstream eout;
	      stringstream ss;
	      ss << "badSequences/badSeq_" << (i + 1) << ".error";
	      efname.assign(ss.str());
	      eout.open(efname.c_str());
	      eout << nameList.at(i) << endl;
	      for (unsigned j = 0; j < seqList.at(i).size(); j++) {
		eout << seqList.at(i).at(j);
	      }
	      eout << endl;
	      eout << "Number of Errors found in this sequence: "
		   << errorCountList.at(i) << endl;
	      eout << errorList.at(i);
	      eout.close();
	    }
	  }
	  fout.close();
	  countGood = nameList.size() - errorCount;
	  cout << "The number of good sequences: " << countGood << endl;
	  cout << "The number of bad sequences: " << errorCount;
	  if (errorCount > 0)
	    cout << ", see the directory 'badSequences'";
	  cout << endl << endl;
	}
	return countGood;
}

