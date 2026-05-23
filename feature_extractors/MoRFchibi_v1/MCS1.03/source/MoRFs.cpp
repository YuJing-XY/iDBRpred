/*
 * MoRFs.cpp
 *
 *  Created on: 2015-09-30
 *      Author: nmalhis
 */

#include "MoRFs.h"
#include <iostream>
#include <fstream>
#include <ostream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <string>

MoRFs::MoRFs() {
	// TODO Auto-generated constructor stub
}

MoRFs::~MoRFs() {
	// TODO Auto-generated destructor stub
}
void MoRFs::loadAnnotation(unsigned int df) {
	unsigned int number_of_sequences;
	if (df < z_files.size()) {
		string line;
		ifstream zin;
		unsigned int p;
		int itmp;
		char ctmp;
		if (z_files.at(df).size() < 2) {
			cerr << "MoRFs::loadAnnotation fatal error: no file is assigned"
					<< endl;
			exit(EXIT_FAILURE);
		}
		zin.open(z_files.at(df).c_str());
		if (!zin.is_open()) {
			cerr << "MoRFs::loadAnnotation fatal error: can\'t open file: "
					<< z_files.at(df) << endl;
			exit(EXIT_FAILURE);
		}
		zin >> number_of_sequences;
		if (number_of_sequences != mdata.at(df).size()) {
			cerr
					<< "MoRFs::loadAnnotation fatal error: not the right annotation file "
					<< z_files.at(df) << "\t" << mdata.at(df).size() << " != "
					<< number_of_sequences << endl;
			exit(EXIT_FAILURE);
		}
		getline(zin, line);
		getline(zin, line);
		while (!zin.eof()) {
			stringstream ss(line);
			ss >> p;
			if (p >= mdata.at(df).size()) {
				cerr << "MoRFs::loadAnnotation fatal error: bad p: " << p
						<< " must be < " << mdata.at(df).size() << endl;
				exit(EXIT_FAILURE);
			}
			ss >> ctmp;
			mdata.at(df).at(p).cl.push_back(ctmp);
			ss >> itmp;
			mdata.at(df).at(p).mSize.push_back(itmp);
			ss >> itmp;
			mdata.at(df).at(p).lfStart.push_back(itmp);
			ss >> itmp;
			mdata.at(df).at(p).mStart.push_back(itmp);
			ss >> itmp;
			mdata.at(df).at(p).rfEnd.push_back(itmp);
			getline(zin, line);
		}
		zin.close();
		labled.at(df) = true;
	} else {
		cerr << "MoRFs::loadAnnotation fatal error: bad df: " << df << endl;
	}
}

void MoRFs::load(unsigned int df) {
	if (df < f_files.size()) {
		string line;
		ifstream fin;
		if (f_files.at(df).size() < 2) {
			cerr << "load fatal error: no file is assigned" << endl;
			exit(EXIT_FAILURE);
		}
		fin.open(f_files.at(df).c_str());
		if (!fin.is_open()) {
			cerr << "load fatal error: can\'t open file: " << f_files.at(df)
					<< endl;
			exit(EXIT_FAILURE);
		}
		MProtein mp;
		getline(fin, line);
		while (!fin.eof()) {
			if (line.at(0) == '>') {
				if (mp.name.size() > 1) {
					mdata.at(df).push_back(mp);
					mp.pSeq.clear();
					mp.name.clear();
				}
				mp.name.assign(line);
			} else {
				for (unsigned int i = 0; i < line.size(); i++) {
					if (line.at(i) != ' ' && line.at(i) != 13)
						mp.pSeq.push_back(line.at(i));
				}
			}
			getline(fin, line);
		}
		if (mp.name.size() > 1) {
			mdata.at(df).push_back(mp);
		}
		fin.close();
	} else {
		cerr << "Load fatal error: bad df: " << df << endl;
	}
}

unsigned int MoRFs::addFasta(string datafile, bool lbld) {

	string zfile;
	string ffile;
	ffile.assign(datafile);
	zfile.assign(datafile.substr(0, datafile.size() - 6));
	zfile.append(".annotation");
	labled.push_back(lbld);
	f_files.push_back(ffile);
	z_files.push_back(zfile);

	vector<MProtein> tmorf;
	mdata.push_back(tmorf);
	return labled.size() - 1;
}

