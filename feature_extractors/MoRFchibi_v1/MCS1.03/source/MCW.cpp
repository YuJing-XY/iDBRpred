/*
 * MC.cpp
 *
 *  Created on: 2015-09-30
 *      Author: nmalhis
 */

#include "MCW.h"

#include <fstream>
#include "svm.h"
#include <stdlib.h>

struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;
struct svm_node *x;             // Nawar

MCW::MCW() {
	x = (struct svm_node *) malloc(10 * sizeof(struct svm_node));
	max_nr_attr = 640;
	cross_validation = nr_fold = 0;
	Local_path.assign("");
	ESpritz_path.assign("");
	PSIBLAST_path.assign("");
	db[0].assign("");
	db[1].assign("");
	mType = _AV;
	C = 500;
	G = 1;
	kType = _RBF;
	PRB = 1;
	mFile.assign("");
	sFile.assign("");
	rFile.assign("");
	fFile.assign("");
	trainingFile.assign("");
	pssm_data[0].assign("");
	pssm_data[1].assign("");
}

void MCW::setPathes(string &local, string &espritz, string &psiblast,
		string &spdb, string &urdb) {
	Local_path.assign(local);
	ESpritz_path.assign(espritz);
	PSIBLAST_path.assign(psiblast);
	db[0].assign(spdb);
	db[1].assign(urdb);
	string aaiFile;
	aaiFile.assign(Local_path);
	aaiFile.append(PATH_DATA);
	aaiFile.append("AAI1");
	aai.load(aaiFile);
	mFile.assign(Local_path);
	sFile.assign(Local_path);
	rFile.assign(Local_path);
	fFile.assign(Local_path);
	trainingFile.assign(Local_path);
	pssm_data[0].assign(Local_path);
	pssm_data[1].assign(Local_path);

	mFile.append(PATH_SVM);
	sFile.append(PATH_SVM);
	rFile.append(PATH_MC);
	fFile.append(PATH_DATA);
	trainingFile.append(PATH_TMP);
	pssm_data[0].append("PSSM/inDir_sp/");
	pssm_data[1].append("PSSM/inDir_ur/");
}

MCW::~MCW() {
	// TODO Auto-generated destructor stub
}

void MCW::set_parameters(unsigned mt, float c, float g, unsigned kt, int prb) {
	mType = mt;
	C = c;
	G = g;
	kType = kt;
	PRB = prb;
	if (mType == _AV) {
		mFile.append("tModel.svm");
		sFile.append("tScale.txt");
		rFile.append("tResults.txt");
		trainingFile.append("tTraining.txt");
		fFile.append("tFeatures.txt");
	} else {
		mFile.append("sModel.svm");
		sFile.append("sScale.txt");
		rFile.append("sResults.txt");
		trainingFile.append("sTraining.txt");
		fFile.append("sFeatures.txt");
	}
	aai.loadSelected_features_index(fFile);
	max_nr_attr = aai.selectedFeaturesIndex.size();
	x = (struct svm_node *) realloc(x,
			(max_nr_attr + 1) * sizeof(struct svm_node));
}

void MCW::score(string tsFile) {
	if (mType == _AV) {
		cout << "Running MoRFchibi: SVMt" << endl;
	} else {
		cout << "Running MoRFchibi: SVMs" << endl;
	}
	cout.flush();
	loadFasta(tsFile);
	aai.loadScale(sFile);
	model = svm_load_model(mFile.c_str()); // must load scale
	scoreTestData(0);
	free_SVM_Memory();
}

unsigned int MCW::loadFasta(string fl) {
	unsigned int df;
	string sTmp;
	sTmp.assign(Local_path);
	sTmp.append(fl);
	df = morf.addFasta(sTmp, false);
	morf.load(df);
	return df;
}

void MCW::scoreTestData(unsigned df) { // sTyp: _AV/_MX
	param.kernel_type = kType; //RBF;
	param.probability = PRB;
	param.gamma = G;
	param.C = C;

	param.svm_type = C_SVC;
	param.degree = 3;
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	if (df < morf.mdata.size()) {
		aai.loadScale(sFile);
		resultsOut.open(rFile.c_str());
		for (unsigned p = 0; p < morf.mdata.at(df).size(); p++) {
			if (p % 10 == 0) {
				cout << '.';
				cout.flush();
			}
			if (morf.mdata.at(df).at(p).pSeq.size() > 25) {
				extract_test_sequences(df, p);
				aai.colorSelected();
				if (mType == _AV) {
					aai.scale();
				} else {
					aai.scaleLocal();
				}
				scoreColoredSequence();
				saveSequenceResult(df, p);
			} else {
				scoreShort(df, p);
			}
		}
		cout << endl;
		resultsOut.close();
	}
}

void MCW::scoreColoredSequence() {
	double pe;
	for (unsigned i = 0; i < aai.ftrs_selected.size(); i++) {
		unsigned str = aai.ftrs_selected.at(i).str;
		unsigned siz = aai.ftrs_selected.at(i).siz;
		pe = predictOneWindow(aai.ftrs_selected.at(i));
		if (mType == _AV) {
			for (unsigned j = str; j < str + siz; j++) {
				scs.at(j) += pe;
				iscs.at(j)++;}
			}
		else {
			for (unsigned j = str; j < str + siz; j++) {
				if (scs.at(j) < pe)
					scs.at(j) = pe;
			}
		}
	}
	if (mType == _AV) {
		for (unsigned j = 0; j < scs.size(); j++) {
			if (iscs.at(j) > 0)
				scs.at(j) = scs.at(j) / iscs.at(j);
		}
	}
}

void MCW::saveSequenceResult(unsigned df, unsigned p) {
	resultsOut << morf.mdata.at(df).at(p).name << endl;
	for (unsigned i = 0; i < morf.mdata.at(df).at(p).pSeq.size(); i++) {
		resultsOut << i + 1 << '\t' << (char) morf.mdata.at(df).at(p).pSeq.at(i)
				<< '\t' << scs.at(i) << endl;
	}
}

void MCW::scoreShort(unsigned df, unsigned p) {
	resultsOut << morf.mdata.at(df).at(p).name << endl;
	for (unsigned i = 0; i < morf.mdata.at(df).at(p).pSeq.size(); i++) {
		resultsOut << i + 1 << '\t' << (char) morf.mdata.at(df).at(p).pSeq.at(i)
				<< "\t0.5" << endl;
	}
}

double MCW::predictOneWindow(Features &features) {
	int i;
	double pe[5];
	for (i = 0; i < max_nr_attr; i++) {
		x[i].index = features.fIndex.at(i);
		x[i].value = features.fValue.at(i);
	}
	x[i].index = -1;
	// -----------------------------------------------------------------------
	svm_predict_probability(model, x, pe);
	// -----------------------------------------------------------------------
	return pe[0];
}

void MCW::extract_test_sequences(unsigned int df, unsigned p) {
	aai.aatdrs.clear();
	unsigned sizMax = 25;
	AATMoRF aatdr;
	aatdr.p = p;
	aatdr.cl = -1;
	for (unsigned siz = 6; siz < sizMax; siz++) {
		for (unsigned str = 0; str <= morf.mdata.at(df).at(p).pSeq.size() - siz;
				str++) {
			aatdr.m.clear();
			aatdr.f.clear();
			aatdr.str = str;
			aatdr.siz = siz;
			for (unsigned i = str; i < str + siz; i++) {
				aatdr.m.push_back(morf.mdata.at(df).at(p).pSeq.at(i));
			}
			unsigned lfs = 0;
			unsigned rfe = str + siz + F_SIZE;
			if (str > lfs + F_SIZE) {
				lfs = str - F_SIZE;
			}
			if (rfe > morf.mdata.at(df).at(p).pSeq.size()) {
				rfe = morf.mdata.at(df).at(p).pSeq.size();
			}

			for (unsigned i = lfs; i < str; i++) {
				aatdr.f.push_back(morf.mdata.at(df).at(p).pSeq.at(i));
			}
			for (unsigned i = str + siz; i < rfe; i++) {
				aatdr.f.push_back(morf.mdata.at(df).at(p).pSeq.at(i));
			}
			aai.aatdrs.push_back(aatdr);
		}
	}
	scs.clear();
	iscs.clear();
	scs.reserve(morf.mdata.at(df).at(p).pSeq.size());
	iscs.reserve(morf.mdata.at(df).at(p).pSeq.size());
	for (unsigned i = 0; i < morf.mdata.at(df).at(p).pSeq.size(); i++) {
		scs.push_back(0);
		iscs.push_back(0);
	}
}

// --------------------------------------------------------- run

int MCW::runEspritz(string iFile, int dbg) {
	cout << "Running ESpritz_D" << endl;
	cout.flush();
	string cmd;
	cmd.assign("rm ");
	cmd.append(Local_path);
	cmd.append("idp/inDir/* 2>");
	cmd.append(Local_path);
	cmd.append("tmp/stderr.txt");
	int ret = system(cmd.c_str());

	int countFiles = 0;
	string line;
	string seq;
	string pName;
	string pfile;
	ifstream fin;
	string sTmp;
	sTmp.assign(Local_path);
	sTmp.append(iFile);
	fin.open(sTmp.c_str());
	ofstream fout0;
	string plist;
	string inf;
	plist.assign(Local_path);
	plist.assign("idp/inDir/pList.txt");
	fout0.open(plist.c_str());
	inf.assign(Local_path);
	inf.append("idp/inDir/in");
	seq.clear();
	getline(fin, line);
	while (!fin.eof()) {
		if (line.size() > 0) {
			if (line.at(0) == '>') {
				if (seq.size() > 10) {
					stringstream ss(stringstream::in | stringstream::out);
					ss << inf;
					ss << countFiles;
					ss << ".fasta";
					pfile.assign(ss.str());
					ofstream fout;
					fout.open(pfile.c_str());
					fout << pName << endl;
					fout0 << pName << endl;
					fout << seq << endl;
					fout.close();
					seq.clear();
					countFiles++;
				}
				pName.assign(line);
			} else {
				seq.append(line);
			}
		}
		getline(fin, line);
	}
	if (seq.size() > 10) {
		stringstream ss(stringstream::in | stringstream::out);
		ss << inf;
		ss << countFiles;
		ss << ".fasta";
		pfile.assign(ss.str());
		ofstream fout;
		fout.open(pfile.c_str());
		fout << pName << endl;
		fout0 << pName << endl;
		fout << seq << endl;
		fout.close();
		countFiles++;
	}
	fout0.close();
	ret = chdir(ESpritz_path.c_str());
	if (dbg)
		cerr << "<---- Debug message> ESpritz directory: " << ESpritz_path
				<< endl;
	string jnk;
	jnk.assign(Local_path);
	jnk.append("tmp/d_junk");
	cmd.assign("./espritz.pl ");
	cmd.append(Local_path);
	cmd.append("idp/inDir D 0 > ");
	cmd.append(jnk);
	cmd.append(" 2> ");
	cmd.append(jnk);
	if (dbg)
		cerr << "<---- Debug message> ESpritz command: " << cmd << endl;
	ret = system(cmd.c_str());
	ret = chdir(Local_path.c_str());
	return 0;
}

int MCW::runPSIBLAST(string iFile, int dbt, int noth, int dbg) {
	if (dbt == _SPDB || dbt == _URDB) {
		cout << "Running PSIBLAST: Aligning to " << db[dbt] << endl;
		cout.flush();
		string path_pList;
		string fls;
		string cmd;
		ofstream fout0;

		cmd.assign("rm ");
		cmd.append(pssm_data[dbt]);
		cmd.append("* 2>");
		cmd.append(Local_path);
		cmd.append("tmp/stderr.txt");
		int ret = system(cmd.c_str());
		cmd.clear();

		path_pList.assign(pssm_data[dbt]);
		path_pList.append("pList.txt");

		fout0.open(path_pList.c_str());
		int countFiles = 0;
		string line;
		string seq;
		string pName;
		string pfile;
		ifstream fin;
		fin.open(iFile.c_str());
		seq.clear();
		getline(fin, line);
		while (!fin.eof()) {
			if (line.size() > 0) {
				if (line.at(0) == '>') {
					if (seq.size() > 10) {
						stringstream ss(stringstream::in | stringstream::out);
						ss.str("");
						ss << pssm_data[dbt];
						ss << "in";
						ss << countFiles;
						ss << ".fasta";
						pfile.assign(ss.str());
						ofstream fout;
						fout.open(pfile.c_str());
						fout << pName << endl;
						fout0 << pName << endl;
						fout << seq << endl;
						fout.close();
						seq.clear();
						countFiles++;
					}
					pName.assign(line);
				} else {
					seq.append(line);
				}
			}
			getline(fin, line);
		}
		if (seq.size() > 10) {
			stringstream ss(stringstream::in | stringstream::out);
			ss.str("");
			ss << pssm_data[dbt];
			ss << "in";
			ss << countFiles;
			ss << ".fasta";
			pfile.assign(ss.str());
			ofstream fout;
			fout.open(pfile.c_str());
			fout << pName << endl;
			fout0 << pName << endl;
			fout << seq << endl;
			fout.close();
			countFiles++;
		}
		fout0.close();

		for (int i = 0; i < countFiles; i++) {
			if (i % 10 == 0) {
				cout << '.';
				cout.flush();
			}
			stringstream ss3(stringstream::in | stringstream::out);
			ss3.str("");
			ss3 << PSIBLAST_path;
			ss3 << "psiblast -query ";
			ss3 << pssm_data[dbt];
			ss3 << "in";
			ss3 << i;
			ss3 << ".fasta -db " << db[dbt];
			ss3 << " -out tmp/out -out_ascii_pssm ";
			ss3 << pssm_data[dbt];
			ss3 << "in";
			ss3 << i;
			if (dbt == _SPDB) {
				ss3
						<< ".pssm -num_iterations 2 -comp_based_stats 0 -num_threads ";
			} else {
				ss3 << ".pssm -num_iterations 2 -num_threads ";
			}
			ss3 << noth;
			ss3 << " 2> ";
			ss3 << Local_path;
			ss3 << "tmp/stderr.txt";
			cmd = ss3.str();
			if (dbg)
				cerr << "<---- Debug message> PSIBLAST command: " << cmd
						<< endl;
			ret = system(cmd.c_str());
		}
		cout << endl;
		cout.flush();
	}
	return 0;
}

void MCW::free_SVM_Memory() {
	svm_free_and_destroy_model(&model);
	svm_destroy_param(&param);
	free(prob.y);
	free(prob.x);
	free(x_space);
}

