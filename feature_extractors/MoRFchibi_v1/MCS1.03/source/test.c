/*
 * test.cpp
 *
 *  Created on: 2015-09-30
 *      Author: nmalhis
 */
#include "AAIndex.h"
#include "MoRFs.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "defs.h"
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "Properties.h"

using namespace std;

void test(Properties &prop, char arg) {
  string cmd;
  int ret;
  if (arg == 'e') {
    cerr << "\nTesting ESpritz" << endl;
    cmd.assign("rm ");
    cmd.append(prop.Local_path);
    cmd.append("idp/inDir/* 2>");
    cmd.append(prop.Local_path);
    cmd.append("tmp/stderr.txt");   
    ret = system(cmd.c_str());
    
    cmd.assign("cp ");
    cmd.append(prop.Local_path);
    cmd.append("inData/p53.fasta ");
    cmd.append(prop.Local_path);
    cmd.append("idp/inDir/");
    ret = system(cmd.c_str());
    
    ret = chdir( prop.ESpritz_path.c_str());
    string jnk;
    jnk.assign(prop.Local_path);
    jnk.append("tmp/d_junk");
    cmd.assign("./espritz.pl ");
    cmd.append(prop.Local_path);
    cmd.append("idp/inDir D 0 > ");
    cmd.append(jnk);
    cmd.append(" 2> ");
    cmd.append(jnk);
    ret = system(cmd.c_str());
    ret = chdir(prop.Local_path.c_str());
    ret = system("ls -l idp/inDir/");
  } else {
    string pdb;
    string pssmDir;
    if (arg == 'u') {
      cerr << "\nTesting PSIBLAST with UniRef90" << endl;
      pssmDir.assign("PSSM/inDir_ur/");
      pdb.assign(prop.UniRef90);
    } else {
      cerr << "\nTesting PSIBLAST with SwissProt" << endl;
      pssmDir.assign("PSSM/inDir_sp/");
      pdb.assign(prop.SwissProt);
    }
    cmd.assign("rm ");
    cmd.append(prop.Local_path);
    cmd.append(pssmDir);
    cmd.append("* 2>");
    cmd.append(prop.Local_path);
    cmd.append("tmp/stderr.txt");   
    ret = system(cmd.c_str());
    
    cmd.assign("cp ");
    cmd.append(prop.Local_path);
    cmd.append("inData/p53.fasta ");
    cmd.append(prop.Local_path);
    cmd.append(pssmDir);
    ret = system(cmd.c_str());
    
    stringstream ss3(stringstream::in | stringstream::out);
    ss3.str("");
    ss3 << prop.PSIBLAST_path;
    ss3 << "psiblast -query ";
    ss3 << prop.Local_path;
    ss3 << pssmDir;
    ss3 << "p53";
    ss3 << ".fasta -db " << pdb;
    ss3 << " -out tmp/out -out_ascii_pssm ";
    ss3 << prop.Local_path;
    ss3 << pssmDir;
    ss3 << "p53.pssm -num_iterations 2 -comp_based_stats 2 -num_threads ";
    ss3 << prop.threads;
    ss3 << " 2> ";
    ss3 << prop.Local_path;
    ss3 << "tmp/stderr.txt";
    cmd.assign(ss3.str());
    ret = system(cmd.c_str());
    cmd.assign("ls -l ");
    cmd.append(pssmDir);
    ret = system(cmd.c_str());
  }
}

int main(int argc, char** argv) {
  int ret;
  ret = system("clear");

  Properties prop;
  prop.soft.assign("Testing the MoRFchibi SYSTEM supporting software installation");
  // prop.makeLocalPath();
  // prop.removeLocalPath(false);
  if (argc == 2) { 
    if (!(argv[1][1] == 0 && (argv[1][0] == 'e' || argv[1][0] == 's' || 
			      argv[1][0] == 'u'))) {
      cerr << "\nSYNOPSIS\n\t./test [e|s|u]\n" << endl;
      cerr << "\nDESCRIPTION\n\t'e' for testing ESpritz installation," << endl;
      cerr << "\t's' for testing PSIBLAST with SwisProt," << endl;
      cerr << "\t'u' for testing PSIBLAST with SwisProt, or" << endl;
      cerr << "\tno argument for all.\n" << endl;
      exit(0);
    }
    char arg = argv[1][0];
    test(prop, arg);
  } else if (argc == 1) {
    char argList[] = {'e', 's', 'u'};
    for (int i = 0; i < 3; i++) {
      test(prop, argList[i]);
    }
  } else {
    cerr << "\nSYNOPSIS\n\t./test [e|s|u]" << endl;
    cerr << "\nDESCRIPTION\n\t'e' for testing ESpritz installation," << endl;
    cerr << "\t's' for testing PSIBLAST with SwisProt," << endl;
    cerr << "\t'u' for testing PSIBLAST with SwisProt, or" << endl;
    cerr << "\tno argument for all.\n" << endl;
    exit(0);
  }
}
