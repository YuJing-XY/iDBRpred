/*
 * Properties.h
 *
 *  Created on: Nov 13, 2015
 *      Author: nmalhis
 */

#ifndef PROPERTIES_H_
#define PROPERTIES_H_

#include <string>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <errno.h>
#include <stdio.h>

using namespace std;

struct Properties {
	int threads;
	int debug;
	string SwissProt;
	string UniRef90;
	string inputFile;
	string outputFile;
	string release;
	string soft;
	string ESpritz_path;
	string PSIBLAST_path;
	string Local_path;


	Properties();

	int clean_fasta(string &cleanOut, unsigned minSize);

	void message(bool disp) {
		cout << "\n\n---\tThe University of British Columbia" << endl;
		cout
				<< "---\tMichael Smith Laboratories - Center for High-Throughput Biology"
				<< endl;
		cout << "---\tGsponer Lab" << endl;
		cout << "" << endl;
		cout << "---\t" << soft << " " << release << endl;
		cout << "---\tDeveloped by: Nawar Malhis" << endl;
		cout << endl;
		cout << "Input file: " << inputFile << endl;
		cout << "Output file: " << outputFile << endl << endl;
		cout.flush();
		if (disp) {
			cout << Local_path << endl;
			cout << ESpritz_path << endl;
			cout << PSIBLAST_path << endl;
			cout << inputFile << endl;
			cout << outputFile << endl;
			cout << SwissProt << endl;
			cout << UniRef90 << endl;
			cout << threads << endl << endl;
		}
	}
};

#endif /* PROPERTIES_H_ */
