#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>


#include "Highs.h"

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
  double objective_function_value, rsumpinf, rsumdinf, sinf, mx_inf, highs_q;

  double *dobj, *dclo, *dcup, *drlo, *drup, *dels, *dels2, *s;
  double *d, *lam, *clp, *lamrow;
  int *sv_bs_r, *sv_bs_c;
  int *colpts, *rix, *mrow, *mrow2, *mcol, *collen;
  char buffer[100];
  string line;

  bool set_stat;

  ifstream inFile;

  HighsBasis basis;
  HighsStatus return_status;
  HighsModelStatus lp_status, lp_status_scaled;

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
  
  printf("About to create Highs instance\n"); fflush(stdout);

  Highs highs;
  HighsLp lp;
  printf("Created HighsLp instance\n"); fflush(stdout);
  lp.numCol_ = n;
  lp.numRow_ = m;
  lp.colCost_.resize(n);
  lp.colLower_.resize(n);
  lp.colUpper_.resize(n);
  lp.Astart_.resize(n+1);
  lp.rowLower_.resize(m);
  lp.rowUpper_.resize(m);
  lp.Aindex_.resize(nel);
  lp.Avalue_.resize(nel);
  
  printf("About to copy col data\n"); fflush(stdout);
  
  for(i=0;i<n;i++) {
    lp.colCost_[i] = dobj[i];
    lp.colLower_[i] = dclo[i];
    lp.colUpper_[i] = dcup[i];
    lp.Astart_[i] = colpts[i];
  }
  printf("Copied col data\n"); fflush(stdout);
  lp.Astart_[n] = colpts[n];
  printf("Copied colpts[n]\n"); fflush(stdout);
  for(i=0;i<m;i++) {
    lp.rowLower_[i] = drlo[i];
    lp.rowUpper_[i] = drup[i];
  }
  printf("Copied row data\n"); fflush(stdout);
  for(i=0;i<nel;i++){
    lp.Aindex_[i] = mrow2[i];
    lp.Avalue_[i] = dels2[i];
  }
  printf("Copied matrix data\n"); fflush(stdout);
  //  flag = highs.loadModel(n, m, nel, dobj, dclo, dcup, drlo, drup, colpts, mrow2, dels2);
  return_status = highs.passModel(lp);

  //-------------------- set warmstart information ----------------
  if (set_stat) {
    //    basis = highs.getBasis();
    basis.col_status.resize(n);
    for (int i = 0; i < basis.col_status.size(); i++) {
      basis.col_status[i] = (HighsBasisStatus)sv_bs_c[i];
    }
    
    basis.row_status.resize(m);
    for (int i = 0; i < basis.row_status.size(); i++) {
      basis.row_status[i] = (HighsBasisStatus)sv_bs_r[i];
    }
    highs.setBasis(basis);
  }

  // Set HiGHS options: in src/lp_data/HighsOptions.h
  highs.setHighsOptionValue("primal_feasibility_tolerance", RTOLPINF);
  highs.setHighsOptionValue("simplex_iteration_limit", IMAXITER);
  if (set_stat)
    highs.setHighsOptionValue("presolve", "off");
  //  highs.setHighsOptionValue("simplex_scale_strategy", 3);

    //  highs.writeModel("ml.mps");
  return_status = highs.run();
  lp_status = highs.getModelStatus(false); //original (unscaled) lp
  lp_status_scaled = highs.getModelStatus(true);  //scaled lp
  HighsInfo info = highs.getHighsInfo();

  printf("Unscaled status: %s\n", highs.highsModelStatusToString(lp_status).c_str());
  printf("Scaled status  : %s\n", highs.highsModelStatusToString(lp_status_scaled).c_str());

  //highs.setHighsOptionValue("message_level", 0);

  const HighsInfo& highs_info = highs.getHighsInfo();
  objective_function_value = highs_info.objective_function_value;
  nr_lp = highs_info.simplex_iteration_count;

  const HighsSolution& solution = highs.getSolution();
  if (solution.col_value.size()) {
    for(i=0;i<n;i++){
      d[i]= solution.col_value[i];
      lam[i]= solution.col_dual[i];
    }
    for(i=0;i<m;i++){
      clp[i] = solution.row_value[i];
      lamrow[i] = solution.row_dual[i];
    }
  //    highs.getSolution(d, lam, clp, lamrow);
   
    for(i=0;i<n;i++){
      if (d[i]-dcup[i]>1e-5  || dclo[i]-d[i] > 1e-5) {
	//      cout << i << " " << dclo[i] << " " << d[i] << " " << dcup[i] << "\n";
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
  
    printf("HiGHS gives Num Primal Infeasibilities = %d\n", highs_info. num_primal_infeasibilities);
    printf("HiGHS gives Max Primal Infeasibility   = %g\n", highs_info. max_primal_infeasibility);
    printf("HiGHS gives Sum Primal Infeasibilities = %g\n", highs_info. sum_primal_infeasibilities);
    printf("HiGHS gives Num   Dual Infeasibilities = %d\n", highs_info. num_dual_infeasibilities);
    printf("HiGHS gives Max   Dual Infeasibility   = %g\n", highs_info. max_dual_infeasibility);
    printf("HiGHS gives Sum   Dual Infeasibilities = %g\n", highs_info. sum_dual_infeasibilities);

    printf("#primal inf = %d\n", inumpinf);
    printf("#    suminf = %f\n", rsumpinf);
    printf("#    maxinf = %f\n", mx_inf);
  }
  printf("#   lp_iter = %d\n", nr_lp);

}


