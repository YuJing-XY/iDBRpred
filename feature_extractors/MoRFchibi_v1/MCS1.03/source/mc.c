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
  prop.soft.assign("MoRFchibi_Light");
  cleanInputFile.assign("clean_");
  cleanInputFile.append(prop.inputFile);
  prop.message(false);
  good_sequences = prop.clean_fasta(cleanInputFile, 15);
  if (good_sequences > 0) {
    Stat st;
    st.soft.assign(prop.soft);
    st.setPathes(prop.Local_path, prop.ESpritz_path, prop.PSIBLAST_path);
    // prop.inputFile.assign(st.validate(prop.inputFile));
    //*
    {
      MCW mcw;
      mcw.setPathes(prop.Local_path, prop.ESpritz_path, prop.PSIBLAST_path, prop.SwissProt, prop.UniRef90);
      mcw.set_parameters(_AV, 500, 1, _RBF, 1);
      mcw.score(cleanInputFile);
    }
    
    {
      MCW mcw;
      mcw.setPathes(prop.Local_path, prop.ESpritz_path, prop.PSIBLAST_path, prop.SwissProt, prop.UniRef90);
      mcw.set_parameters(_MX, 500, 0.001, _Sigmoid, 1);
      mcw.score(cleanInputFile);
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
    
    st.addEmpty(6);
    
    // ------------- 8 MC
    st.load_Joined2(0, 1, "MC");
    st.reScale_load();
    
    
    vector<unsigned> vroc;
    vroc.push_back(8);
    vroc.push_back(1);
    vroc.push_back(0);
    st.output(vroc, prop.outputFile, prop.release);
  }
  //*/
}
