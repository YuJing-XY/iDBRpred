/*
 * MoRFs.h
 *
 *  Created on: 2015-09-30
 *      Author: nmalhis
 */

#ifndef MORFS_H_
#define MORFS_H_

#include <string>
#include <vector>
#include <stdlib.h>
#include <iostream>
// #include <pthread.h>
// #include <chrono>
#include <unistd.h>
#include "defs.h"

using namespace std;

struct MProtein {
	string name;
	vector<int> pSeq;
	vector<int> label;

	vector<char> cl;
	vector<int> mSize;
	vector<int> mStart;
	vector<int> lfStart;
	vector<int> rfEnd;
};

class MoRFs {
public:
	vector<vector<MProtein> > mdata;
	vector<bool> labled;
	vector<string> m_files;
	vector<string> f_files;
	vector<string> z_files;

	MoRFs();
	virtual ~MoRFs();

	void loadAnnotation(unsigned int df);
	void load(unsigned int df);
	unsigned int addFasta(string datafile, bool lbld);
};

#endif /* MORFS_H_ */
