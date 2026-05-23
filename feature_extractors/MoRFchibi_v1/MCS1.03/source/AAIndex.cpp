/*
 * AAIndex.cpp
 *
 *  Created on: 2015-09-30
 *      Author: nmalhis
 */

#include "AAIndex.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
// #include "svm.h"

using namespace std;

AAIndex::AAIndex() {
	// TODO Auto-generated constructor stub
	int taa[] = { 0, 17, 13, 3, 2, 16, 4, 6, 7, 8, 11, 10, 12, 5, 15, 18, 19,
			22, 24, 21 };
	for (int i = 0; i < 20; i++) {
		i_to_aa[i] = taa[i];
	}
	int tpIndex[] = { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5,
			1, 15, 16, -1, 19, 17, -1, 18 };
	for (int i = 0; i < 25; i++) {
		aa_to_i[i] = tpIndex[i];
	}
	saveScale = true;
}

AAIndex::~AAIndex() {
	// TODO Auto-generated destructor stub
}

int AAIndex::load(const string& sourceFile) {
	AAIRec oneRec;
	char flag;
	string line;
	string word;
	ifstream fin;
	this->sourceFile = sourceFile;
	if (sourceFile.size() < 2) {
		cerr << "Fatal error: no AAIndex file is assigned" << endl;
		exit(EXIT_FAILURE);
	}
	fin.open(sourceFile.c_str());
	if (!fin.is_open()) {
		cerr << "Fatal error: can\'t open AAIndex file: " << sourceFile << endl;
		exit(EXIT_FAILURE);
	}

	getline(fin, line);
	while (!fin.eof()) {
		stringstream ss(line);
		ss >> word;
		if (word.size() < 3 && word.size() > 0) {
			flag = word.at(0);
		}
		if (flag == '/' && oneRec.vals.size() == 20) {
			aai1_data.push_back(oneRec);
			oneRec.vals.clear();
			flag = ' ';
		} else if (flag == 'H') {
			ss >> oneRec.accessionNumber;
			flag = ' ';
		} else if (flag == 'D') {
			line.erase(0, 2);
			oneRec.dataDescription = line;
			flag = ' ';
		} else if (flag == 'C') {
		} else if (flag == 'I') {
			float tmp;
			float max = -10000;
			float avrg = 0;
			for (int j = 0; j < 2; j++) {
				getline(fin, line);
				stringstream ss(line);
				for (int i = 0; i < 10; i++) {
					ss >> tmp;
					avrg += tmp;
					oneRec.vals.push_back(tmp);
				}
			}
			avrg = avrg / 20;
			for (unsigned int i = 0; i < oneRec.vals.size(); i++) {
				oneRec.vals.at(i) = oneRec.vals.at(i) - avrg;
				if (fabs(oneRec.vals.at(i)) > max)
					max = fabs(oneRec.vals.at(i));
			}
			for (unsigned int i = 0; i < oneRec.vals.size(); i++) {
				oneRec.vals.at(i) = oneRec.vals.at(i) / max;
			}
			flag = ' ';
		} else {
			flag = ' ';
		}
		line.clear();
		getline(fin, line);
	}
	fin.close();
	return 0;
}

void AAIndex::loadSelected_features_index(string fFile) {
	string line;
	ifstream fin;
	int iTmp;
	selectedFeaturesIndex.clear();
	fin.open(fFile.c_str());
	if (!fin.is_open()) {
		cerr << "load fatal error: can\'t open file: " << fFile << endl;
		exit(EXIT_FAILURE);
	}
	fin >> iTmp;
	while (!fin.eof()) {
		if (iTmp < 0 || iTmp >= ((int) aai1_data.size() * 2)) {
			cerr << "Error: bad feature index: " << iTmp << " in file " << fFile
					<< " at index " << selectedFeaturesIndex.size() << endl;
		} else {
			selectedFeaturesIndex.push_back(iTmp);
		}
		fin >> iTmp;
	}
}

void AAIndex::colorSelected() {
	int fidx;
	float av;
	int ai;
	ftrs_selected.clear();
	for (unsigned int i = 0; i < aatdrs.size(); i++) { // for each training sample
		Features sfs; // selected features
		sfs.fIndex.clear();
		sfs.fValue.clear();
		sfs.cv = -1;
		sfs.p = aatdrs.at(i).p;
		sfs.cl = aatdrs.at(i).cl;
		sfs.str = aatdrs.at(i).str;
		sfs.siz = aatdrs.at(i).siz;
		for (unsigned int fi = 0; fi < selectedFeaturesIndex.size(); fi++) {
			av = 0;
			fidx = selectedFeaturesIndex.at(fi) / 2;
			if (selectedFeaturesIndex.at(fi) % 2 == 0) {
				for (unsigned int j = 0; j < aatdrs.at(i).m.size(); j++) {
					ai = aa_to_i[aatdrs.at(i).m.at(j) - 'A'];
					av += aai1_data.at(fidx).vals.at(ai);
				}
				av = av / aatdrs.at(i).m.size();
			} else {
				for (unsigned int j = 0; j < aatdrs.at(i).f.size(); j++) {
					ai = aa_to_i[aatdrs.at(i).f.at(j) - 'A'];
					av += aai1_data.at(fidx).vals.at(ai);
				}
				av = av / aatdrs.at(i).f.size();
			}
			sfs.fValue.push_back(av);
			sfs.fIndex.push_back(selectedFeaturesIndex.at(fi));
		}
		ftrs_selected.push_back(sfs);
	}
}

void AAIndex::genScale(string oFile) {
	float av;
	float mx;
	_av.clear();
	_mx0.clear();
	for (unsigned i = 0; i < selectedFeaturesIndex.size(); i++) {
		av = 0;
		mx = -1000;
		for (unsigned j = 0; j < ftrs_selected.size(); j++) {
			av += ftrs_selected.at(j).fValue.at(i);
		}
		av = av / ftrs_selected.size();

		for (unsigned j = 0; j < ftrs_selected.size(); j++) {
			if (fabs(ftrs_selected.at(j).fValue.at(i) - av) > mx) {
				mx = fabs(ftrs_selected.at(j).fValue.at(i) - av);
			}
		}
		if (mx == 0)
			cerr << "AAIndex::scale. Error: MX must not be zero: "
					<< _mx0.size() << endl;
		_av.push_back(av);
		_mx0.push_back(mx);
	}
	if (saveScale) {
		ofstream fout;
		fout.open(oFile.c_str());
		for (unsigned i = 0; i < selectedFeaturesIndex.size(); i++) {
			fout << _av.at(i) << '\t' << _mx0.at(i) << endl;
		}
		fout.close();
	}
}

void AAIndex::loadScale(string iFile) {
	ifstream fin;
	float av;
	float mx;
	string line;
	_av.clear();
	_mx0.clear();
	fin.open(iFile.c_str());
	for (unsigned i = 0; i < selectedFeaturesIndex.size(); i++) {
		getline(fin, line);
		if (!fin.eof()) {
			stringstream ss(line);
			ss >> av;
			ss >> mx;
			if (mx == 0)
				cerr << "AAIndex::scale. Error (load) : MX must not be zero: "
						<< _mx0.size() << endl;
			_av.push_back(av);
			_mx0.push_back(mx);
		} else {
			cerr << "AAIndex::scale. Error: EOF\t" << iFile << endl;
		}
	}
	fin.close();
}

void AAIndex::scale() {
	if (_mx0.size() == selectedFeaturesIndex.size()
			&& _av.size() == selectedFeaturesIndex.size()) {
		for (unsigned i = 0; i < selectedFeaturesIndex.size(); i++) {
			for (unsigned j = 0; j < ftrs_selected.size(); j++) {
				ftrs_selected.at(j).fValue.at(i) =
						(ftrs_selected.at(j).fValue.at(i) - _av.at(i))
								/ _mx0.at(i);
			}
		}
	} else {
		cerr << "AAIndex::scale() Error";
	}
}

void AAIndex::scaleLocal() {
	float av;
	float mx0;
	for (unsigned i = 0; i < selectedFeaturesIndex.size(); i++) {
		av = 0;
		mx0 = -1000;
		for (unsigned j = 0; j < ftrs_selected.size(); j++) {
			av += ftrs_selected.at(j).fValue.at(i);
		}
		av = av / ftrs_selected.size();

		for (unsigned j = 0; j < ftrs_selected.size(); j++) {
			ftrs_selected.at(j).fValue.at(i) -= av;
			if (fabs(ftrs_selected.at(j).fValue.at(i)) > mx0)
				mx0 = fabs(ftrs_selected.at(j).fValue.at(i));
		}
		if (mx0 == 0)
			cerr << "AAIndex::scale. Error: MX must not be zero: "
					<< i << endl;
		for (unsigned j = 0; j < ftrs_selected.size(); j++) {
			ftrs_selected.at(j).fValue.at(i) = ftrs_selected.at(j).fValue.at(i) / mx0;
		}
	}
}

void AAIndex::saveTrainingData(string trFile, int cvOut) {
	ofstream fout;
	fout.open(trFile.c_str());
	if (!fout.is_open()) {
		cerr << "AAIndex::saveTrainingData fatal error: can\'t open file: "
				<< trFile << endl;
		exit(EXIT_FAILURE);
	}
	for (unsigned int i = 0; i < ftrs_selected.size(); i++) {
		if (ftrs_selected.at(i).cv != cvOut) {
			fout << (char) ftrs_selected.at(i).cl;
			for (unsigned int j = 0; j < ftrs_selected.at(i).fIndex.size();
					j++) {
				fout << "\t" << ftrs_selected.at(i).fIndex.at(j) << ":"
						<< ftrs_selected.at(i).fValue.at(j);
			}
			fout << endl;
		}
	}
	fout.close();
}
