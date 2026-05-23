/*
 * trun.cpp
 *
 *  Created on: 2015-09-30
 *      Author: nmalhis
 */

#include <iostream>

#include "MCW.h"
#include "Stat.h"
#include "Properties.h"

using namespace std;

int main(int argc, char** argv) {
  string cleanInputFile;
  int good_sequences = 0;
  Properties prop;
  prop.soft.assign("MoRFchibi SYSTEM");
  cleanInputFile.assign("clean_");
  cleanInputFile.append(prop.inputFile);
  prop.message(false);
  good_sequences = prop.clean_fasta(cleanInputFile, 15);
  if (good_sequences > 0) {
    // cout << "The number of good sequences: " << good_sequences << endl;
    Stat st;
    st.soft.assign(prop.soft);
    st.setPathes(prop.Local_path, prop.ESpritz_path, prop.PSIBLAST_path);
    // prop.inputFile.assign(st.validate(prop.inputFile));
    
    {
      MCW mcw;
      mcw.setPathes(prop.Local_path, prop.ESpritz_path,
		    prop.PSIBLAST_path, prop.SwissProt, prop.UniRef90);
      mcw.set_parameters(_AV, 500, 1, _RBF, 1);
      mcw.score(cleanInputFile);
    }
    
    {
      MCW mcw;
      mcw.setPathes(prop.Local_path, prop.ESpritz_path,
		    prop.PSIBLAST_path, prop.SwissProt, prop.UniRef90);
      mcw.set_parameters(_MX, 500, 0.001, _Sigmoid, 1);
      mcw.score(cleanInputFile);
    }

    {
      MCW mcw;
      mcw.setPathes(prop.Local_path, prop.ESpritz_path,
		    prop.PSIBLAST_path, prop.SwissProt, prop.UniRef90);
      mcw.runEspritz(cleanInputFile, prop.debug);
    }
    
    {
      MCW mcw;
      mcw.setPathes(prop.Local_path, prop.ESpritz_path,
		    prop.PSIBLAST_path, prop.SwissProt, prop.UniRef90);
      mcw.runPSIBLAST(cleanInputFile, _SPDB, prop.threads, prop.debug);
    }
    
    {
      MCW mcw;
      mcw.setPathes(prop.Local_path, prop.ESpritz_path,
		    prop.PSIBLAST_path, prop.SwissProt, prop.UniRef90);
      mcw.runPSIBLAST(cleanInputFile, _URDB, prop.threads, prop.debug);
    }
    
    string svmtResults;
    string svmsResults;
    svmtResults.assign(prop.Local_path);
    svmtResults.append(PATH_MC);
    svmsResults.assign(svmtResults);
    svmtResults.append("tResults.txt");
    svmsResults.append("sResults.txt");
    
    st.loadScores_MC(svmtResults, "SVMt"); // -------------- 0 SVMt
    st.loadScores_MC(svmsResults, "SVMs"); // -------------- 1 SVMs
    st.load_Espritz(); // ------------- 2 ESpritz_D
    st.reScale_load();
    st.load_PSSM(1);
    st.reScale_load();
    st.load_Reverse(); // ------------- 4 RWGMtP_Reverse
    st.load_PSSM(2);   // ------------- 5 Info, 6 WOP
    st.reScale_load();
    st.reScale_load(-2);
    // ------------- 7 ICS
    set<int> ics_comps;
    ics_comps.insert(-1);
    ics_comps.insert(-2);
    ics_comps.insert(-3);
    st.load_Joined(ics_comps, "ICS");
    st.reScale_load();
    
    // ------------- 8 MC
    st.load_Joined2(0, 1, "MC");
    st.reScale_load();
    
    // ------------- 9 MCS
    st.load_UpDown(7, 2);
    st.reScale_load();
    
    // ------------- 10 MDC
    st.load_Joined2(9, 2, "MDC");
    st.reScale_load();
    
    // ------------- 11 MCW
    st.load_Joined2(10, 8, "MCW");
    
    // ------------- 12 MCL
    st.load_Joined2(2, 8, "MCL");
    
    vector<unsigned> vroc;
    vroc.push_back(11);
    vroc.push_back(12);
    vroc.push_back(8);
    vroc.push_back(10);
    vroc.push_back(2);
    vroc.push_back(7);
    st.output(vroc, prop.outputFile, prop.release);
  } else {
    cerr << "Error: No good sequences in " << prop.inputFile 
	 << ", see the directory 'badSequences'" << endl;
  }
}
