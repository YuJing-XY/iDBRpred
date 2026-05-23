/*
 * Stat.h
 *
 *  Created on: Oct 30, 2015
 *      Author: nmalhis
 */

#ifndef STAT_H_
#define STAT_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>
#include <set>
#include <cmath>
#include <algorithm>
#include "defs.h"

using namespace std;

class PRef {
public:
	int idx;
	int first;
	int afterLast;
	string name;
};

class VRef {
public:
	set<PRef>::iterator it;
	int cluster;
};

class Comp {
public:
	bool operator()(PRef p1, PRef p2) const {
		if (p1.name.compare(p2.name) < 0)
			return true;
		else
			return false;
	}
};

class Stat {
 public:
  set<PRef, Comp> pId;
  vector<set<PRef>::iterator> seqList;
  float _auc[3];
  float _wauc[3];
  vector<vector<float> > pe;
  vector<int> cl;
  vector<char> AA;
  int upe;
  float distribution[2][1001];
  float fptp[2][1001];
  float roc[1001];
  vector<string> dataSource;
  string Local_path;
  string ESpritz_path;
  string PSIBLAST_path;
  string pssm_data[2];
  string soft;

  // MCW vars
  string rsName;
  int upRadius;
  int downRadius;
  float upPsThreshold;
  float downPsThreshold;
  int pCountThreshold;
  float coverUpCut;
  float coverDownCut;
  int sqCount[7];
  vector<float> rDistribution;
  vector<int> Up;
  vector<int> Down;
  vector<float> pMax;
  vector<float> AVCover;

  // -------------------------------------
  Stat();
  virtual ~Stat();
  string validate(string &ifile);
  void addEmpty(int x);
  float bayes(float f1, float f2);
  void reScale_load(int d = -1);                    // MCW ok
  void normalize1(unsigned &d, float &mx, float &mn);   // MCW ok
  void fillGaussian(float sigRG = 1);               // MCW ok
  void load_Joined2(int d1, int d2, string nm);
  void load_Joined(set<int> &comps, string nm);     // MCW
  void loadScores_MC(string scFile, string dsName);
  void load_Espritz();                              // MCW
  void load_PSSM(int cnt);                          // MCW
  void load_Reverse(int source = -1);               // MCW
  void load_UpDown(int source, int cover);
  void output(vector<unsigned> &vpe, string oFile, string release);
  void loadAnnotation(string anFile);
  void saveAnnotation(string junkFile);
  float compRoc(int dpe, unsigned clSz);
  void compRoc(vector<unsigned> &vpe);
  void mark_Up(int source, int cover);
  void mark_Down(int source, int cover);
  void markLongAnnotation(int lngStart);
  void setPathes(string &local, string &espritz, string &psiblast);

  void resetUpDown() {                             // MCW
    Up.clear();
    Down.clear();
    pMax.clear();
    AVCover.clear();
    for (unsigned i = 0; i < AA.size(); i++) {
      Up.push_back(0);
      Down.push_back(0);
      pMax.push_back(0);
      AVCover.push_back(0);
    }
  }

  void listProteins() {
    for (set<PRef, Comp>::iterator it = pId.begin(); it != pId.end();
	 ++it) {
      cerr << it->name << "\t" << it->first << "\t" << it->afterLast
	   << endl;
    }
    for (unsigned i = 0; i < seqList.size(); i++) {
      cerr << seqList.at(i)->name << "\t" << seqList.at(i)->first << "\t"
	   << seqList.at(i)->afterLast << endl;
    }
  }

  int getAAindex(char aa) {                         // MCW
    switch (aa) {
    case 'A':
      return 0;
    case 'R':
      return 1;
    case 'N':
      return 2;
    case 'D':
      return 3;
    case 'C':
      return 4;
    case 'Q':
      return 5;
    case 'E':
      return 6;
    case 'G':
      return 7;
    case 'H':
      return 8;
    case 'I':
      return 9;
    case 'L':
      return 10;
    case 'K':
      return 11;
    case 'M':
      return 12;
    case 'F':
      return 13;
    case 'P':
      return 14;
    case 'S':
      return 15;
    case 'T':
      return 16;
    case 'W':
      return 17;
    case 'Y':
      return 18;
    case 'V':
      return 19;
    }
    cerr << "Error Stat:getAAindex. Bad aa: " << aa << endl;
    return -1;
  }

  void countUpDown() {
    int countUp = 0;
    int countDown = 0;
    for (unsigned i = 0; i < AA.size(); i++) {
      if (Up.at(i) > 0) {
	countUp++;
      }
      if (Down.at(i) > 0) {
	countDown++;
      }
    }
    cerr << countUp << "\t" << countDown << endl;
  }

  unsigned get_ed(int dd) {
    unsigned ed = abs(dd);
    if (ed >= pe.size()) {
      cerr << "Error in Stat::get_ed: bad ed: " << ed << "\n ed is set to -1" << endl;
      ed = pe.size() - 1;
    } else if (dd < 0) {
      ed = pe.size() - ed;
    }
    return ed;
  }
};

#endif /* STAT_H_ */
