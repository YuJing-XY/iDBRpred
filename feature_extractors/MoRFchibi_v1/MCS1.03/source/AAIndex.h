/*
 * AAIndex.h
 *
 *  Created on: 2015-09-30
 *      Author: Nawar Malhis
 */

#ifndef AAINDEX_H_
#define AAINDEX_H_

#include <string>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include "defs.h"

using namespace std;

struct AAIRec {
	vector<float> vals;
	string accessionNumber;
	string dataDescription;
};

struct AATMoRF { // amino acid training data record
	vector<char> m;
	vector<char> f;
	int cl;
	int p;
	unsigned str;
	unsigned siz;
};

struct Features {
	int cv;
	int p;
	int cl;
	unsigned str;
	unsigned siz;
	vector<int> fIndex;
	vector<float> fValue;
};

class AAIndex {
public:
	string path_trFile;
	string sourceFile;
	bool saveScale;
	vector <float> _av;
	vector <float> _mx0;
	vector<int> selectedFeaturesIndex;
	vector<AAIRec> aai1_data;
	vector<AATMoRF> aatdrs;
	vector <Features> ftrs_selected;

	int i_to_aa[20]; // to convert AA into int
	// = {0, 17, 13, 3, 2, 16, 4, 6, 7, 8, 11, 10, 12, 5, 15, 18, 19, 22, 24, 21};
	int aa_to_i[26]; // to convert int into AA (aa in reverse)
	// = {0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18};

	AAIndex();
	virtual ~AAIndex();

	int load(const string& sourceFile);
	void loadSelected_features_index(string sFile);
	void colorSelected();
	void genScale(string oFile);
	void loadScale(string iFile);
	void scale();
	void scaleLocal();
	void saveTrainingData(string trFile, int cvOut);

 	void displayRec(unsigned int idx) {
		if (idx >= 0 && idx < aai1_data.size()) {
			cout << aai1_data.at(idx).accessionNumber << endl;
			cout << aai1_data.at(idx).dataDescription << endl;
			for (unsigned int i = 0; i < aai1_data.at(idx).vals.size(); i++) {
				cout << (char) ('A' + i_to_aa[i]) << "\t"
						<< aai1_data.at(idx).vals.at(i) << endl;
			}
		}
	}

 	void displayAVMX() {
		for (unsigned int i = 0; i < selectedFeaturesIndex.size(); i++) {
			cout << i << "\t" << _av.at(i) << "\t" << _mx0.at(i) << endl;
		}
	}

 	void displaySelectedFeaturesIndex() {
		for (unsigned int i = 0; i < selectedFeaturesIndex.size(); i++) {
			cout << i << "\t" << selectedFeaturesIndex.at(i) << endl;
		}
	}
};

#endif /* AAINDEX_H_ */
