#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>


#include "highs_c_api.h"
#include "HConst.h"
#include "Highs.h"
#include "HighsStatus.h"
#include "HighsLp.h"

// ------------------------ some option values --------------------------

#define RTOLPINF 1e-5
//#define RTOLPINF 1e-7
//#define RTOLPINF 9.9999999999999991e-06

#define IMAXITER 10000

using namespace std;


int
main(int argc, char *argv[])
{
  int n, m, nel;
  int i;
  int nr_lp, inumpinf, inumdinf, iprobstat, flag, highs_stat;
  double q, rsumpinf, rsumdinf, sinf, mx_inf, highs_q;

  double *dobj, *dclo, *dcup, *drlo, *drup, *dels, *dels2, *s;
  double *d, *lam, *clp, *lamrow;
  int *sv_bs_r, *sv_bs_c;
  int *colpts, *rix, *mrow, *mrow2, *mcol, *collen;
  char buffer[100];
  string line;

  bool set_stat;

  ifstream inFile;

  void *highs;
  HighsBasis basis;

  // ------------------------ read data file  --------------------------
  // (in matlab format with hexadecimal floating point number format)

  set_stat = false;

  printf("Open File: %s\n", argv[1]);
  //inFile.open("lp0036.m");
  inFile.open(argv[1]);

  getline(inFile, line);
  getline(inFile, line);
  strncpy(buffer, line.c_str(), 99);
  n = atoi(buffer);
  getline(inFile, line);
  getline(inFile, line);
  getline(inFile, line);
  strncpy(buffer, line.c_str(), 99);
  m = atoi(buffer);
  getline(inFile, line);
  getline(inFile, line);
  getline(inFile, line);
  strncpy(buffer, line.c_str(), 99);
  nel = atoi(buffer);
  getline(inFile, line);
  getline(inFile, line);
  
  dclo = new double[n];
  dcup = new double[n];
  dobj = new double[n];
  collen = new int[n];
  colpts = new int[n+1];
  drlo = new double[m];
  drup = new double[m];
  s = new double[n+m];
  sv_bs_c = new int[n];
  sv_bs_r = new int[m];
  dels = new double[nel];
  dels2 = new double[nel];
  mrow = new int[nel];
  mrow2 = new int[nel];
  mcol = new int[nel];

  d = new double[n];
  lam = new double[n];
  lamrow = new double[m];
  clp = new double[m];



  for (i=0;i<n;i++){
    //getline(inFile, line);
    //strncpy(buffer, line.c_str(), 99);
    //strtok(buffer, ",; \n\r");
    //dclo[i] = atof(strtok(NULL, ",; \n\r"));
    //dcup[i] = atof(strtok(NULL, ",; \n\r"));
    //dobj[i] = atof(strtok(NULL, ",; \n\r"));

    getline(inFile, line);
    strncpy(buffer, line.c_str(), 99);
    strtok(buffer, ",; \n\r");

    sscanf(strtok(NULL, ",; \n\r"), "%la", &(dclo[i]));
    sscanf(strtok(NULL, ",; \n\r"), "%la", &(dcup[i]));
    sscanf(strtok(NULL, ",; \n\r"), "%la", &(dobj[i]));

    //printf("%d %f %f %f\n", i+1, dclo[i], dcup[i], dobj[i]);
  }
  getline(inFile, line);
  getline(inFile, line);

  for (i=0;i<m;i++){
    getline(inFile, line);
    strncpy(buffer, line.c_str(), 99);
    strtok(buffer, ",; \n\r");
    //drlo[i] = atof(strtok(NULL, ",; \n\r"));
    //drup[i] = atof(strtok(NULL, ",; \n\r"));
    sscanf(strtok(NULL, ",; \n\r"), "%la", &(drlo[i]));
    sscanf(strtok(NULL, ",; \n\r"), "%la", &(drup[i]));
  }
  getline(inFile, line);
  getline(inFile, line);

  for (i=0;i<nel;i++){
    getline(inFile, line);
    strncpy(buffer, line.c_str(), 99);
    strtok(buffer, ",; \n\r");
    mrow[i] = atoi(strtok(NULL, ",; \n\r"));
    mcol[i] = atoi(strtok(NULL, ",; \n\r"));
    mrow[i]--;
    mcol[i]--;

    //dels[i] = atof(strtok(NULL, ",; \n\r"));
    sscanf(strtok(NULL, ",; \n\r"), "%la", &(dels[i]));
  }
  getline(inFile, line);
  getline(inFile, line);

  getline(inFile, line);
  if (!inFile.eof()){
    strncpy(buffer, line.c_str(), 99);
    highs_stat = atoi(strtok(buffer, ",; \n\r"));
    highs_q = atof(strtok(NULL, ",; \n\r"));
  }else{
    highs_stat = 0;
    highs_q = 0;
  }
  inFile.close();

  set_stat = false;
  if (argc>2){
    printf("Open File: %s\n", argv[2]);
    inFile.open(argv[2]);

    getline(inFile, line);
    for(i=0;i<n;i++){
      getline(inFile, line);
      sv_bs_c[i] = atoi(line.c_str());
    }     
    for(i=0;i<m;i++){
      getline(inFile, line);
      sv_bs_r[i] = atoi(line.c_str());
    }     

    inFile.close();
    set_stat = true;
  }

 //  Need to convert row wise order to column wise order
  
  // count number of elements in each column
  for(i=0;i<n;i++) collen[i] = 0;
  for (i=0;i<nel;i++) collen[mcol[i]]++;
  colpts[0] = 0;
  for(i=0;i<n-1;i++) colpts[i+1] = colpts[i] + collen[i];
  // and go fill in all entries
  for(i=0;i<nel;i++){
    dels2[colpts[mcol[i]]] = dels[i];
    mrow2[colpts[mcol[i]]] = mrow[i];
    colpts[mcol[i]]++;
  }
  colpts[0] = 0;
  for(i=0;i<n;i++) colpts[i+1] = colpts[i] + collen[i];

  // ======================== set up Highs  ========================
  

  highs = Highs_create();
  flag = Highs_loadModel(highs, n, m, nel, dobj, dclo, dcup, drlo, drup,
			   colpts, mrow2, dels2);
 


  //-------------------- set warmstart information ----------------
  if (set_stat) {
    basis = ((Highs *)highs)->getBasis();
    for (int i = 0; i < basis.col_status.size(); i++) {
      basis.col_status[i] = (HighsBasisStatus)sv_bs_c[i];
    }
    
    for (int i = 0; i < basis.row_status.size(); i++) {
      basis.row_status[i] = (HighsBasisStatus)sv_bs_r[i];
    }
    ((Highs *)highs)->setBasis(basis);
    //Highs_setBasis(highs, basis);
  }      
  

  // Set HiGHS options: in src/lp_data/HighsOptions.h
  sprintf(buffer, "%g", RTOLPINF);
  Highs_setHighsOptionValue(highs, "primal_feasibility_tolerance", buffer);
  sprintf(buffer, "%d", IMAXITER);
  Highs_setHighsOptionValue(highs, "simplex_iteration_limit", buffer);

  iprobstat = Highs_run(highs);
   
  //sprintf(outline, "%d", 0);
  //Highs_setHighsOptionValue(highs, "message_level", outline);

  Highs_getSolution(highs, d, lam, clp, lamrow);
   
  q = Highs_getObjectiveValue(highs);
  nr_lp = Highs_getIterationCount(highs);

  Highs_destroy(highs);

  for(i=0;i<n;i++){
    if (d[i]-dcup[i]>1e-5  || dclo[i]-d[i] > 1e-5) {
      cout << i << " " << dclo[i] << " " << d[i] << " " << dcup[i] << "\n";
    }
    //cout << i << " " << lam(i) << "\n";
    //normd = max(normd,fabs(d[i]/s[i]));
  }


 for(i=0;i<n;i++){
    // FIXME: s[i] was dspace(ncolscales-1+i). Is that correct?
    mx_inf = max(mx_inf, max(
			     (d[i]-dcup[i])*s[i],
			     (dclo[i]-d[i])*s[i]));
    //cout << i << " " << s[i] << "\n";
  }
  for (i=0;i<m;i++){
    // FIXME: s[n+i] was dspace(nrowscales-1+i). Is that correct?
    mx_inf = max(mx_inf, max(
			     (clp[i]-drup[i])*s[n+i],
			     (drlo[i]-clp[i])*s[n+i]));
    //cout << i << " " << s[n+i] << "\n";
  }
  
  inumpinf = 0;
  rsumpinf = 0.0;
  for(i=0;i<n;i++){
    sinf =  max(max(d[i]-dcup[i], dclo[i]-d[i]),0.);
    //sinf = scalar infeasibility
    if (sinf>RTOLPINF) inumpinf++;
    rsumpinf += sinf;
  }
  for(i=0;i<m;i++){
    sinf =  max(max(clp[i]-drup[i], drlo[i]-clp[i]),0.);
    if (sinf>RTOLPINF) inumpinf++;
    rsumpinf += sinf;
  }
  
  printf("#primal inf = %d\n", inumpinf);
  printf("#    suminf = %f\n", rsumpinf);
  printf("#    maxinf = %f\n", mx_inf);

}


