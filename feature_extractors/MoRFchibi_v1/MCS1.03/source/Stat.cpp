/*
 * Stat.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: nmalhis
 */

#include "Stat.h"

Stat::Stat() {
	rsName.assign("");
	upe = -1;
	upRadius = 3;           // 3
	downRadius = 7;         // 7
	pCountThreshold = 3;    // 3
	upPsThreshold = 0.45;   // 0.45
	downPsThreshold = 0.60; // 0.60
	coverUpCut = 0.45;      // 0.45
	coverDownCut = 0.45;    // 0.45
	sqCount[0] = 0;
	resetUpDown();
	for (int i = 1; i < 7; i++)
		sqCount[i] = 2;
	pssm_data[0].assign("");
	pssm_data[1].assign("");
}

void Stat::setPathes(string &local, string &espritz, string &psiblast) {
	Local_path.assign(local);
	ESpritz_path.assign(espritz);
	PSIBLAST_path.assign(psiblast);
	rsName.assign(Local_path);
	rsName.append(PATH_MC);
	rsName.append("reScale/rescale_0.txt");
	pssm_data[0].assign(Local_path);
	pssm_data[1].assign(Local_path);
	pssm_data[0].append("PSSM/inDir_sp/");
	pssm_data[1].append("PSSM/inDir_ur/");
}

Stat::~Stat() {
	// TODO Auto-generated destructor stub
	int ret = system("rm tmp/*");
}

string Stat::validate(string &ifile) {
	ifstream fin;
	string outFile;
	string ch;
	ch.assign("A");
	outFile.assign("tmp/input.fasta");
	fin.open(ifile.c_str());
	if (!fin.is_open()) {
		cerr << "Stat::validate fatal error: can\'t open file: " << ifile
				<< endl;
		exit(EXIT_FAILURE);
	}
	string line;
	int linenumber;
	string outData;
	// char aa;
	int ret;
	linenumber = 1;
	outData.clear();
	getline(fin, line);
	while (!fin.eof()) {
		if (outData.size() == 0) {
			if (line.size() != 0)
				if (line.at(0) != '>') {
					cerr << "Stat::validate fatal error: bad fasta, line:\n\t"
							<< line << endl;
					cerr << " must start with '>'" << endl;
					exit(0);
				}
		}
		if (line.at(0) == '>') {
			if (outData.size() > 1)
				outData.append("\n");
			outData.append(line);
			outData.append("\n");
		} else {
			if (line.size() > 0) {
				for (unsigned i = 0; i < line.size() - 1; i++) {
					ret = getAAindex(line.at(i));
					if (ret >= 0) {
						// ch.assign("A");
						ch.at(0) = line.at(i);
						outData.append(ch);
					} else {
						cerr << "Stat::validate fatal error: bad amino acid ("
								<< line.at(i) << ")" << endl;
						cerr << "\tin line: " << linenumber << " at index: "
								<< (i + 1) << endl;
						cerr << line << endl;
						exit(0);
					}
				}
				ret = getAAindex(line.at(line.size() - 1));
				if (ret >= 0) {
					ch.at(0) = line.at(line.size() - 1);
					outData.append(ch);
				}
			}
		}
		linenumber++;
		getline(fin, line);
	}
	if (outData.size() > 1)
		outData.append("\n");
	fin.close();
	ofstream fout;
	fout.open(outFile.c_str());
	for (unsigned i = 0; i < outData.size(); i++) {
		fout << outData.at(i);
	}
	fout.close();
	return outFile;
}

void Stat::addEmpty(int x) {
	vector<float> vf;
	vf.clear();
	for (int i = 0; i < x; i++) {
		pe.push_back(vf);
		dataSource.push_back("Not Used");
	}
}

float Stat::bayes(float f1, float f2) {
	return ((f1 * f2) / ((f1 * f2) + ((1 - f1) * (1 - f2))));
}

void Stat::reScale_load(int dd) {
	unsigned ed = get_ed(dd);
	float mx;
	float mn;
	float ftmp;
	vector<float> map;
	map.clear();
	ifstream fin;
	rsName.at(rsName.size() - 5) = char('0' + ed);
	fin.open(rsName.c_str());
	fin >> mx;
	fin >> mn;
	for (int j = 0; j < DistributionSIZE; j++) {
		fin >> ftmp;
		map.push_back(ftmp);
	}
	fin.close();

	normalize1(ed, mx, mn);
	for (unsigned i = 0; i < pe.at(ed).size(); i++) {
		pe.at(ed).at(i) = map.at(
				int(pe.at(ed).at(i) * (DistributionSIZE - 1) + 0.49));
		if (pe.at(ed).at(i) > 0.95)
			pe.at(ed).at(i) = 0.95;
		if (pe.at(ed).at(i) < 0.05)
			pe.at(ed).at(i) = 0.05;
	}
}

void Stat::normalize1(unsigned &ed, float &mx, float &mn) {
	if (ed >= 0 && ed < pe.size()) {
		if (mx == -1 && mn == -1) {
			mx = -100;
			mn = 100;
			for (unsigned i = 0; i < pe.at(ed).size(); i++) {
				if (pe.at(ed).at(i) > mx)
					mx = pe.at(ed).at(i);
				if (pe.at(ed).at(i) < mn)
					mn = pe.at(ed).at(i);
			}
			mx = mx - mn;
		}
		for (unsigned i = 0; i < pe.at(ed).size(); i++) {
			pe.at(ed).at(i) = (pe.at(ed).at(i) - mn) / mx;
			pe.at(ed).at(i) = pe.at(ed).at(i) * 0.9 + 0.05;
			if (pe.at(ed).at(i) > 0.996)
				pe.at(ed).at(i) = 0.996;
			if (pe.at(ed).at(i) < 0.004)
				pe.at(ed).at(i) = 0.004;
		}
	}
}

void Stat::fillGaussian(float sigRG) {
	float sig = (DistributionSIZE - 1) / 10;
	float mu = (DistributionSIZE - 1) / 2;
	float p1;
	float gam;
	if (sigRG > 0 && sigRG < 1)
		sig = sig * sigRG;
	p1 = (float) 1.0 / (sig * sqrt(2 * 3.1418));
	rDistribution.clear();
	for (int i = 0; i < DistributionSIZE; i++) {
		gam = p1 * exp(((mu - i) * (i - mu)) / (2 * sig * sig));
		rDistribution.push_back(gam);
	}
	for (int i = 1; i < DistributionSIZE; i++) {
		rDistribution.at(i) = rDistribution.at(i) + rDistribution.at(i - 1);
	}
}

void Stat::load_Joined2(int dd1, int dd2, string nm) {
	string title;
	int d1 = get_ed(dd1);
	int d2 = get_ed(dd2);
	title.assign(nm);
	dataSource.push_back(title);
	int vAll = pe.size();
	vector<float> vf;
	vf.clear();
	pe.push_back(vf);
	for (unsigned i = 0; i < cl.size(); i++) {
		pe.back().push_back(0.5);
	}
	for (unsigned i = 0; i < cl.size(); i++) {
		if (pe.at(d1).at(i) > 0 && pe.at(d1).at(i) < 1 && pe.at(d2).at(i) > 0
				&& pe.at(d2).at(i) < 1) {
			pe.at(vAll).at(i) = bayes(pe.at(d1).at(i), pe.at(d2).at(i));
		} else {
			cerr << "Error in Stat::load_Joined2 input out of range at:" << i
					<< "\t" << pe.at(d1).at(i) << "\t" << pe.at(d2).at(i)
					<< endl;
		}
		if (pe.at(vAll).at(i) > 1 || pe.at(vAll).at(i) <= 0)
			cerr << "Error in Stat::load_Joined2: pe out of range: "
					<< pe.at(vAll).at(i) << endl;
	}
}

void Stat::load_Joined(set<int> &comps, string nm) {
	string title;
	set<unsigned> cps;
	for (set<int>::iterator it = comps.begin(); it != comps.end(); ++it) {
		cps.insert(get_ed(*it));
	}

	title.assign(nm);
	dataSource.push_back(title);
	int vAll = pe.size();
	vector<float> vf;
	vf.clear();
	pe.push_back(vf);
	for (unsigned i = 0; i < AA.size(); i++) {
		pe.back().push_back(0.5);
	}
	for (unsigned i = 0; i < AA.size(); i++) {
		for (int j = 0; j < vAll; j++) {
			if (cps.find(j) != cps.end()) {
				if (pe.at(j).at(i) > 0 && pe.at(j).at(i) < 1) {
					pe.at(vAll).at(i) = bayes(pe.at(vAll).at(i),
							pe.at(j).at(i));
				} else {
					cerr << "Error in Stat::load_Joined: pe.at(" << j << ").at("
							<< i << "): " << pe.at(j).at(i) << endl;
				}
			}
		}
		if (pe.at(vAll).at(i) > 1 || pe.at(vAll).at(i) <= 0)
			cerr << "Error in Stat::load_Joined: pe out of range: "
					<< pe.at(vAll).at(i) << endl;
	}
}

void Stat::loadScores_MC(string mcFile, string dsName) {
	PRef pref;
	int counter;
	set<PRef>::iterator it;
	vector<float> fv;
	ifstream zin;
	zin.open(mcFile.c_str());
	cout << "Loading scores: " << dsName << endl;
	cout.flush();
	if (!zin.is_open()) {
		cerr << "Stat::loadScores_MC fatal error: can\'t open file: " << mcFile
				<< endl;
		exit(EXIT_FAILURE);
	}
	dataSource.push_back(dsName);
	pe.push_back(fv);
	if (pe.size() == 1) {
		AA.clear();
		cl.clear();
	}
	upe = pe.size() - 1;
	string line;
	int idx;
	char aa;
	float sc;
	pref.first = -1;
	counter = 0;
	getline(zin, line);
	while (!zin.eof()) {
		if (line.at(0) == '>') {
			if (pe.size() == 1) {
				if (pref.first >= 0) {
					pref.afterLast = cl.size();
					if (!pId.insert(pref).second) {
						cerr << "Stat::loadScores_MC Error: " << pref.name
								<< " found" << endl;
					} else {
						seqList.push_back(it);
					}
				}
				pref.idx = counter;
				pref.first = cl.size();
				if (line.at(line.size() - 1) < 30) {
					line.assign(line.substr(0, line.size() - 1));
				}
				pref.name.assign(line);
				counter++;
			} else {
				if (line.at(line.size() - 1) < 30) {
					line.assign(line.substr(0, line.size() - 1));
				}
				pref.name.assign(line);
				it = pId.find(pref);
				if (it == pId.end()) {
					cerr << "Stat::loadScores_MC Error: " << pref.name
							<< " not found" << endl;
				}
			}
		} else {
			stringstream ss(line);
			ss >> idx;
			ss >> aa;
			ss >> sc;
			pe.back().push_back(sc);
			if (pe.size() == 1) {
				cl.push_back(0);
				AA.push_back(aa);
			}
		}
		getline(zin, line);
	}
	if (pe.size() == 1) {
		pref.afterLast = cl.size();
		if (!pId.insert(pref).second) {
			cerr << "Stat::loadScores_MC last Error: " << pref.name << " found"
					<< endl;
		} else {
			seqList.push_back(it);
		}
	}
	for (set<PRef, Comp>::iterator it = pId.begin(); it != pId.end(); ++it) {
		seqList.at(it->idx) = it;
	}
}

void Stat::load_Espritz() {
	dataSource.push_back("IDP");
	int countFiles;
	ifstream fin0;
	ifstream fin;
	string tmp;
	string line;
	string pName;

	int location;
	float sc;
	PRef pref;
	set<PRef>::iterator it;
	vector<float> vf;
	pe.push_back(vf);
	cout << "Loading scores: ESpritz_D" << endl;
	cout.flush();
	for (unsigned i = 0; i < AA.size(); i++) {
		pe.back().push_back(0.5);
	}
	string plist;
	string inf;
	plist.assign(Local_path);
	plist.append("idp/inDir/pList.txt");
	inf.assign(Local_path);
	inf.append("idp/inDir/in");
	fin0.open(plist.c_str());
	getline(fin0, line);
	pref.name.assign(line);
	countFiles = 0;
	while (!fin0.eof()) {
		it = pId.find(pref);
		if (it == pId.end()) {
			cerr << "Stat::load_Espritz_All Error: " << pref.name
					<< " not found" << endl;
			location = -1;
		} else {
			location = it->first;
			string fName;
			stringstream ss0(stringstream::in | stringstream::out);
			ss0.str("");
			ss0 << inf;
			ss0 << countFiles;
			ss0 << ".espritz";
			fName.assign(ss0.str());
			fin.open(fName.c_str());
			countFiles++;
			for (int i = 0; i < 9; i++) {
				getline(fin, line);
			}
			while (!fin.eof()) {
				if (location < it->afterLast) {
					stringstream ss(stringstream::in | stringstream::out);
					ss.str(line);
					ss >> tmp;
					ss >> tmp;
					sc = (float) atof(tmp.c_str());
					pe.back().at(location) = sc;
					location++;
				} else {
					cerr << "Stat::load_Espritz_All: Error: after last: "
							<< location << "\t" << it->first << "\t"
							<< it->afterLast << "\t" << it->name << endl;
				}
				getline(fin, line);
			}
			fin.close();
		}
		getline(fin0, line);
		pref.name.assign(line);
	}
}

void Stat::load_PSSM(int cnt) {
	if (cnt == 1 || cnt == 2) {
		int countFiles = 0;
		int pssl[20];
		int wop[20];
		float info;
		float wmtp;
		string pass;
		string tmp;
		string line;
		int LLin;
		char AAin;
		ifstream fin0;
		ifstream fin;
		PRef pref;
		set<PRef>::iterator it;
		vector<float> vf;
		for (unsigned i = 0; i < AA.size(); i++) {
			vf.push_back(0.5);
		}
		if (cnt == 2) {
			dataSource.push_back("Info");
			dataSource.push_back("WOP");
			pe.push_back(vf);
			pe.push_back(vf);
			cout << "Loading scores: PSSM UniRef90" << endl;
		} else if (cnt == 1) {
			dataSource.push_back("RWGMtP");
			pe.push_back(vf);
			cout << "Loading scores: PSSM SwissProt" << endl;
		}
		cout.flush();
		int location = 0;
		string plist;
		plist.assign(pssm_data[cnt - 1]);
		plist.append("pList.txt");
		fin0.open(plist.c_str());
		getline(fin0, line);
		pref.name.assign(line);
		countFiles = 0;
		while (!fin0.eof()) {
			it = pId.find(pref);
			if (it == pId.end()) {
				cerr << pref.name << " not found" << endl;
				location = -1;
			} else {
				location = it->first;
				string fName;
				stringstream ss0(stringstream::in | stringstream::out);
				ss0.str("");
				ss0 << pssm_data[cnt - 1];
				ss0 << "in";
				ss0 << countFiles;
				ss0 << ".pssm";
				fName.assign(ss0.str());
				fin.open(fName.c_str());
				countFiles++;
				for (int i = 0; i < 4; i++) {
					getline(fin, line);
				}
				while (!fin.eof()) {
					if (line.size() > 50) {
						stringstream ss(stringstream::in | stringstream::out);
						ss.str(line);
						ss >> LLin;
						ss >> AAin;
						if (AAin != AA.at(location)) {
							// sleep(1);
							cerr << fName << "\t" << location - it->first
									<< " Bad AA ref: " << AA.at(location)
									<< "\tAAin: " << AAin << endl;
						}
						for (int i = 0; i < 20; i++) {
							ss >> pssl[i]; // position-specific scoring line
						}
						for (int i = 0; i < 20; i++) {
							ss >> wop[i]; // weighted observed percentages
						}
						ss >> info; // information per position
						ss >> wmtp; // weight of matches to pseudocounts
						if (cnt == 2) {
							pe.at(pe.size() - 2).at(location) = info;
							pe.at(pe.size() - 1).at(location) = float(
									wop[getAAindex(AAin)]) / 100;
						} else {
							pe.at(pe.size() - 1).at(location) = wmtp;
						}
						location++;
						if (location > it->afterLast) {
							cerr << "PSSM bad location " << location << " "
									<< pref.name << endl;
						}
					} else {
					}
					getline(fin, line);
				}
				fin.close();
			}
			getline(fin0, line);
			pref.name.assign(line);
		}
		fin0.close();
	}
}

void Stat::output(vector<unsigned>& vpe, string oFile, string release) {
	for (unsigned i = 0; i < vpe.size(); i++) {
		if (vpe.at(i) >= cl.size()) {
			cerr << "Error Stat::output, bad data set index: " << vpe.at(i)
					<< endl;
			return;
		}
	}
	ofstream fout;
	fout.open(oFile.c_str());
	if (!fout.is_open()) {
		cerr << "Error Stat::output: can\'t open file: " << oFile << endl;
	} else {
		fout << "#" << endl;
		fout << "# " << soft << " " << release << endl;
		fout << "#" << endl;
		fout << "# The University of British Columbia" << endl;
		fout
				<< "# Michael Smith Laboratories - Center for High-Throughput Biology "
				<< endl;
		fout << "#" << endl;
		fout << "# Column\tData type" << endl;
		fout << "#   1   \tresidue index" << endl;
		fout << "#   2   \tresidue" << endl;
		for (unsigned j = 0; j < vpe.size(); j++) {
			fout << "#   " << j + 3 << "   \t" << dataSource.at(vpe.at(j))
					<< endl;
		}
		fout << "#" << endl;
		for (unsigned i = 0; i < seqList.size(); i++) {
			fout << seqList.at(i)->name << endl;
			int counter = 1;
			for (int ix = seqList.at(i)->first; ix < seqList.at(i)->afterLast;
					ix++) {
				fout << counter++ << "\t" << AA.at(ix);
				for (unsigned j = 0; j < vpe.size(); j++) {
					fout << "\t" << pe.at(vpe.at(j)).at(ix);
				}
				fout << endl;
			}
		}
	}
	fout.close();
}

void Stat::load_Reverse(int source) {
	if (source < 0)
		source = pe.size() - 1;
	string dSc;
	dSc.assign(dataSource.at(source));
	dSc.append("_Reverse");
	dataSource.push_back(dSc);
	vector<float> vf;
	vf.clear();
	for (unsigned i = 0; i < pe.at(source).size(); i++) {
		if (pe.at(source).at(i) < 1 && pe.at(source).at(i) > 0)
			vf.push_back(1 - pe.at(source).at(i));
		else
			cerr << "Error: " << pe.at(source).at(i) << endl;
	}
	pe.push_back(vf);
}

void Stat::load_UpDown(int source, int cover) {
	float mVal;
	string dSc;
	dSc.assign(dataSource.at(source));
	dSc.append("_UD");
	dataSource.push_back(dSc);
	vector<float> vf;
	vf.clear();
	pe.push_back(vf);
	for (unsigned i = 0; i < AA.size(); i++) {
		pe.back().push_back(pe.at(source).at(i));
	}
	resetUpDown();
	mark_Up(source, cover);
	mark_Down(source, cover);
	for (unsigned i = 0; i < Up.size(); i++) {
		if (Up.at(i) > 0) {
			// ----------------------------------------------
			Down.at(i) = 0;
			mVal = pMax.at(i);
			for (int j = 0; j < Up.at(i); j++) {
				if (Up.at(i) <= 5)
					for (int q = 0; q < sqCount[Up.at(i) - 1]; q++) {
						mVal = sqrt(mVal);
					}
			}
			if (pe.back().at(i) < mVal)
				pe.back().at(i) = mVal;
			// -----------------------------------------------
		} else if (Down.at(i) > 0) {
			pe.back().at(i) = pe.back().at(i) * AVCover.at(i);
		}
	}
}

void Stat::mark_Up(int source, int cover) {
	float avCover;
	float mVal;
	int uCount;
	float pmax;
	unsigned diameter = 1 + 2 * upRadius;
	float uCut = diameter * coverUpCut;

	unsigned i = 0;
	avCover = 0;
	uCount = 0;
	pmax = 0;
	// ---------------------------------------------
	for (; i < diameter; i++) {
		avCover += pe.at(cover).at(i);
		if (pe.at(source).at(i) > upPsThreshold) {
			uCount++;
			pmax += pe.at(source).at(i);
		}
	}

	if (uCount >= pCountThreshold && avCover > uCut) {
		mVal = pmax / uCount;
		for (int j = 0; j <= upRadius; j++) {
			pMax.at(j) = mVal;
			Up.at(j) = uCount - pCountThreshold + 1;
			if (Up.at(j) > 4)
				Up.at(j) = 4;
		}
	}
	// ---------------------------------------------
	for (; i < AA.size(); i++) {
		avCover -= pe.at(cover).at(i - diameter);
		if (pe.at(source).at(i - diameter) > upPsThreshold) {
			uCount--;
			pmax -= pe.at(source).at(i - diameter);
			if (pmax < 0)
				pmax = 0;
		}
		avCover += pe.at(cover).at(i);
		if (pe.at(source).at(i) > upPsThreshold) {
			uCount++;
			pmax += pe.at(source).at(i);
		}

		if (uCount >= pCountThreshold && avCover > uCut) {
			mVal = pmax / uCount;
			pMax.at(i - upRadius) = mVal;
			Up.at(i - upRadius) = uCount - pCountThreshold + 1;
			if (Up.at(i - upRadius) > 4)
				Up.at(i - upRadius) = 4;
		}
		// }
		// ---------------------------------------------
		if (uCount >= pCountThreshold && avCover > uCut) {
			mVal = pmax / uCount;
			for (unsigned j = i + 1 - upRadius; j < i; j++) {
				pMax.at(j) = mVal;
				Up.at(j) = uCount - pCountThreshold + 1;
				if (Up.at(j) > 4)
					Up.at(j) = 4;
			}
		}
	} // ---------------------------------------------
	  //*/
}

void Stat::mark_Down(int source, int cover) {
	float avCover;
	int dCount;
	unsigned diameter = 1 + 2 * downRadius;
	float dCut = diameter * coverDownCut;

	unsigned i = 0;
	avCover = 0;
	dCount = 0;
	// ---------------------------------------------
	for (; i < diameter; i++) {
		avCover += pe.at(cover).at(i);
		if (pe.at(source).at(i) > downPsThreshold) {
			dCount++;
		}
	}
	if (dCount == 0 && avCover < dCut) {
		for (int j = 0; j <= downRadius; j++) {
			Down.at(j) = 1;
			AVCover.at(j) = avCover / diameter;
		}
	}
	// ---------------------------------------------
	for (; i < AA.size(); i++) {
		avCover -= pe.at(cover).at(i - diameter);
		if (pe.at(source).at(i - diameter) > downPsThreshold) {
			dCount--;
		}
		avCover += pe.at(cover).at(i);
		if (pe.at(source).at(i) > downPsThreshold) {
			dCount++;
		}

		if (dCount == 0 && avCover < dCut) {
			Down.at(i - downRadius) = 1;
			AVCover.at(i - downRadius) = avCover / diameter;
		}
	}
	// ---------------------------------------------
	if (dCount == 0 && avCover < dCut) {
		for (unsigned j = i + 1 - downRadius; j < i; j++) {
			Down.at(j) = 1;
			AVCover.at(j) = avCover / diameter;
		}
	}
}

// ---------------------------------------------------- Delete the rest
void Stat::loadAnnotation(string anFile) {
	if (seqList.size() > 0) {
		unsigned int number_of_sequences;
		string line;
		ifstream zin;
		unsigned p;
		unsigned itmp;
		unsigned sztmp;
		unsigned sttmp;
		unsigned ctmp;
		zin.open(anFile.c_str());
		if (!zin.is_open()) {
			cerr << "Stat::loadAnnotation fatal error: can\'t open file: "
					<< anFile << endl;
			exit(EXIT_FAILURE);
		}
		zin >> number_of_sequences;
		if (number_of_sequences != seqList.size()) {
			cerr
					<< "Stat::loadAnnotation fatal error: not the right annotation file "
					<< anFile << " number_of_sequences: " << number_of_sequences
					<< " != " << seqList.size() << endl;
			exit(EXIT_FAILURE);
		}
		getline(zin, line);
		getline(zin, line);
		while (!zin.eof()) {
			stringstream ss(line);
			ss >> p;
			if (p >= pId.size()) {
				cerr << "Stat::loadAnnotation fatal error: bad p: " << p
						<< " must be < " << seqList.size() << endl;
				exit(EXIT_FAILURE);
			}
			ss >> ctmp;
			if (ctmp != 0) {
				ss >> sztmp;
				ss >> itmp;
				ss >> sttmp;
				if (seqList.at(p)->first + sttmp + sztmp <= cl.size()) {
					for (unsigned i = seqList.at(p)->first + sttmp;
							i < seqList.at(p)->first + sttmp + sztmp; i++) {
						cl.at(i) = 1;
					}
				} else {
					cerr << "Stat::annotation outside the cl vector "
							<< seqList.at(p)->first + sttmp + sztmp << " <= "
							<< cl.size() << endl;
					exit(EXIT_FAILURE);
				}
			}
			getline(zin, line);
		}
		zin.close();
	} else {
		cerr << "Stat::loadAnnotation fatal error: bad file: " << anFile
				<< endl;
	}
}

void Stat::saveAnnotation(string junkFile) {
	int x;
	unsigned sIdx = 0;
	ofstream fout;
	fout.open(junkFile.c_str());
	for (unsigned i = 0; i < cl.size(); i++) {
		if (sIdx < seqList.size()) {
			if (i == (unsigned) seqList.at(sIdx)->first) {
				fout << ">name" << endl;
				sIdx++;
			}
		}
		x = i - seqList.at(sIdx - 1)->first;
		fout << x << "\t" << cl.at(i) << endl;
	}
	fout.close();
}

float Stat::compRoc(int dpe, unsigned clSz) {
	float mx;
	float mn;
	int cl_lbl;
	if (clSz > 2)
		clSz = 2;
	unsigned cltg = clSz;
	if (cltg == 0)
		cltg = 3;
	upe = dpe;
	int countCl[2];
	mx = 0;
	mn = 1;
	countCl[0] = countCl[1] = 0;
	for (unsigned i = 0; i < pe.at(upe).size(); i++) {
		if (pe.at(upe).at(i) > mx)
			mx = pe.at(upe).at(i);
		if (pe.at(upe).at(i) < mn)
			mn = pe.at(upe).at(i);
	}

	mx = mx - mn;
	for (unsigned i = 0; i < pe.at(upe).size(); i++) {
		pe.at(upe).at(i) = pe.at(upe).at(i) - mn;
		pe.at(upe).at(i) = pe.at(upe).at(i) / mx;
	}
	for (int i = 0; i < 1001; i++) {
		distribution[0][i] = distribution[1][i] = 0;
	}
	for (unsigned i = 0; i < pe.at(upe).size(); i++) {
		cl_lbl = (cl.at(i) & cltg) > 0;
		if (!(cl_lbl == 0 && cl.at(i) != 0)) {
			distribution[cl_lbl][int(1000 * pe.at(upe).at(i) + 0.5)]++;
			countCl[cl_lbl]++;
		}
	}
	for (int i = 0; i < 1001; i++) {
		distribution[0][i] = distribution[0][i] / countCl[0];
		distribution[1][i] = distribution[1][i] / countCl[1];
	}
	for (int i = 0; i < 1001; i++) {
		fptp[0][i] = fptp[1][i] = 0;
		roc[i] = 0;
	}
	_auc[clSz] = 0;
	_wauc[clSz] = 0;
	fptp[0][1000] = distribution[0][1000];
	fptp[1][1000] = distribution[1][1000];

	for (int i = 999; i >= 0; i--) {
		fptp[0][i] = fptp[0][i + 1] + distribution[0][i];
		fptp[1][i] = fptp[1][i + 1] + distribution[1][i];
	}
	for (int i = 1000; i >= 0; i--) {
		fptp[0][i] = fptp[0][i] / fptp[0][0];
		fptp[1][i] = fptp[1][i] / fptp[1][0];
	}
	int idx;
	int j;
	float step;
	for (int i = 0; i < 1001; i++) {
		idx = int(1000 * fptp[0][1000 - i] + 0.4999);
		if (fptp[1][1000 - i] > 0)
			roc[idx] = fptp[1][1000 - i];
	}
	for (int i = 1; i < 1001; i++) {
		if (roc[i] == 0 && roc[i - 1] > 0) {
			for (j = i + 1; j < 1001; j++) {
				if (roc[j] > 0)
					break;
			}
			step = (roc[j] - roc[i - 1]) / float(j - i + 1);
			for (int k = i; k < j; k++) {
				roc[k] = roc[k - 1] + step;
			}
		}
	}

	for (int i = 0; i < 250; i++) {
		_wauc[clSz] += (((float(250) - float(i)) / float(250)) + 1)
				* (roc[i] / 1001);
		_auc[clSz] += (roc[i] / 1001);
	}

	for (int i = 250; i < 1001; i++) {
		_wauc[clSz] += (roc[i] / 1001);
		_auc[clSz] += (roc[i] / 1001);
	}

	return _wauc[clSz];
}

void Stat::compRoc(vector<unsigned> &vpe) {
	for (unsigned i = 0; i < vpe.size(); i++) {
		if (vpe.at(i) >= cl.size()) {
			cerr << "Error Stat::compRoc, bad data set index: " << vpe.at(i)
					<< endl;
			return;
		}
	}
	for (unsigned i = 0; i < vpe.size(); i++) {
		compRoc(vpe.at(i), 0);
		compRoc(vpe.at(i), 1);
		compRoc(vpe.at(i), 2);
		cout << dataSource.at(vpe.at(i)) << "\t" << _auc[0] << "\t" << _auc[1]
				<< "\t" << _auc[2] << endl;
	}

}

void Stat::markLongAnnotation(int lngStart) {
	int st;
	int sz;
	for (unsigned p = 0; p < seqList.size(); p++) {
		st = seqList.at(p)->first;
		for (int i = seqList.at(p)->first; i < seqList.at(p)->afterLast; i++) {
			if (cl.at(i) == 0) {
				if (i > st) {
					sz = i - st;
					if (sz >= lngStart) {
						for (int j = st; j < i; j++) {
							cl.at(j) = 2;
						}
					}
				}
				st = i + 1;
			}
		}
		if (seqList.at(p)->afterLast > st) {
			sz = seqList.at(p)->afterLast - st;
			if (sz >= lngStart) {
				for (int j = st; j < seqList.at(p)->afterLast; j++) {
					cl.at(j) = 2;
				}
			}
		}
	}
}

