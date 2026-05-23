/*
 * MC.h
 *
 *  Created on: 2015-09-30
 *      Author: nmalhis
 */

#ifndef MCW_H_
#define MCW_H_

#include "AAIndex.h"
#include "MoRFs.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "defs.h"
#include <string.h>
#include <errno.h>
#include <stdio.h>

using namespace std;

class MCW {
public:
	ofstream resultsOut;
	AAIndex aai;
	MoRFs morf;
	int max_nr_attr;                // Nawar
	unsigned mType;
	unsigned kType;
	float C;
	float G;
	int PRB;
	string trainingFile;
	string mFile; // svm model file
	string rFile; // results file
	string sFile; // scale file
	string fFile; // selected features file
	string Local_path;
	string ESpritz_path;
	string PSIBLAST_path;
	string db[2];
	string pssm_data[2];

	vector<int> au;
	vector<int> sau;
	vector<double> scs;
	vector<int> iscs;

	int cross_validation;
	int nr_fold;

	MCW();
	virtual ~MCW();

	void set_parameters(unsigned mt, float c, float g, unsigned kt,
			int prb);
	void score(string tsFile);
	unsigned int loadFasta(string fl);
	void extract_test_sequences(unsigned df, unsigned p);
	void scoreTestData(unsigned df);
	void scoreColoredSequence();
	double predictOneWindow(Features &features);
	void saveSequenceResult(unsigned df, unsigned p);
	void scoreShort(unsigned df, unsigned p);

	int runEspritz(string iFile, int dbg);                      // MCW
	int runPSIBLAST(string iFile, int dbt, int noth, int dbg);  // MCW
	void setPathes(string &local, string &espritz, string &psiblast, string &spdb, string &urdb);

	string inputFasta(string iFile) {
		string path_iFile;
		path_iFile.assign(Local_path);
		path_iFile.append(PATH_DATA);
		path_iFile.append(iFile);
		path_iFile.append(".fasta");
		return path_iFile;
	}
	// ------------------------------------------------------------------------------
	void free_SVM_Memory();
};

#endif /* MCW_H_ */
