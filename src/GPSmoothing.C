#include "WCPLEEANA/GPKernel.h"
#include "WCPLEEANA/GPRegressor.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

//helper class to track indices conversion from 1D to 5D
class CyclicArray {
  public:
    CyclicArray(int nbins[5]) {
      for (int i=0;i<5;i++) {
        n.push_back(nbins[i]);
        x.push_back(0);
      }
    }
    void Increment() { IncrementElement(4); }
    std::vector<int> X() { return x; }
  protected:
    std::vector<int> x,n;
    void IncrementElement(int i) {
      x[i]++;
      if (x[i]==n[i]) {
        x[i] = 0;
        if (i!=0) { IncrementElement(i-1); }
      }
    }
};

//helper function to split a string into a vector of strings, split by specified delimiter
std::vector<std::string> split(std::string line, std::string delimiter){
    std::vector<std::string> list;
    std::string token;
    std::size_t pos = 0;
    while ((pos = line.find(delimiter)) != std::string::npos) {
        token = line.substr(0, pos);
        list.push_back(token);
        line.erase(0, pos + delimiter.length());
    }
    list.push_back(line);
    return list;
}


void GPSmoothing (TVectorD* vec_mean, TMatrixD* cov_mat_bootstrapping, std::string input_filename, bool do_smoothing=true) {

  if (!do_smoothing) { return; }

  int nbins[5];
  std::vector<double> bin_centers_temp, params;
  std::vector<std::vector<double>> bin_centers;

  //read in from config file
  std::string line, temp;
  std::ifstream file;
  file.open(input_filename,std::ios::in);
  if (file.is_open()) {
    //first line specifies the parameters for GPRegresor
    if (!std::getline(file,line)) { return; }
    std::vector v_line = split(line,"\t");
    for (int i=0;i<6;i++) { params.push_back(std::stod(v_line[i])); }
    //next 5 lines specify the bin centers for each dimension
    for (int i=0;i<5;i++) {
      bin_centers_temp.clear();
      if (!std::getline(file,line)) { return; }
      v_line = split(line,"\t");
      nbins[i] = v_line.size();
      for (int j=0;j<nbins[i];j++) { bin_centers_temp.push_back(std::stod(v_line[j])); }
      bin_centers.push_back(bin_centers_temp);
    }
  } else { return; }
  file.close();

  std::vector<GPPoint> gp_points, gp_points_posterior;
  CyclicArray* indices = new CyclicArray(nbins);

  //Create vector of GPPoints
  int nbins_mean = vec_mean->GetNrows();
  for (int i=0;i<nbins_mean;i++) {
    std::vector<int> bin_indices = indices->X();
    indices->Increment();
    double pt_arr[5] = {bin_centers[0][bin_indices[0]], bin_centers[1][bin_indices[1]], bin_centers[2][bin_indices[2]], bin_centers[3][bin_indices[3]], bin_centers[4][bin_indices[4]]};
    gp_points.push_back(GPPoint(pt_arr));
    gp_points_posterior.push_back(GPPoint(pt_arr));
  }

  //Fit and Predict
  RBFKernel kern = RBFKernel(params);
  GPRegressor reg = GPRegressor(kern, true, cov_mat_bootstrapping);		//bool is for normalizing data
  reg.Fit(gp_points, (*vec_mean), false);						//bool is for tuning hyperparameters
  reg.Predict(gp_points_posterior);

  TVectorD vec_mean_temp              = reg.PosteriorMean();
  TMatrixD cov_mat_bootstrapping_temp = reg.PosteriorCov();

  //Debugging -------

  //Rescale cov matrix
  double vec_mean_var = 0;
  double vec_mean_av = 0;
  double vec_mean_max = 0;
  for (int i=0;i<nbins_mean;i++) {
    vec_mean_av += (*vec_mean)[i];
    vec_mean_max = TMath::Max(vec_mean_max,(*vec_mean)[i]);
  }
  vec_mean_av /= nbins_mean;
  for (int i=0;i<nbins_mean;i++) { vec_mean_var += TMath::Power((*vec_mean)[i]-vec_mean_av,2); }
  vec_mean_var /= nbins_mean;
  TMatrixD cov_mat_bootstrapping_scaled(nbins_mean,nbins_mean);
  for (int i=0;i<nbins_mean;i++) {
    for (int j=0;j<nbins_mean;j++) {
      cov_mat_bootstrapping_scaled(i,j) = (*cov_mat_bootstrapping)(i,j)/TMath::Power(vec_mean_var,1);
    }
  }

  TMatrixD correlation_matrix(nbins_mean,nbins_mean);
  TMatrixD correlation_matrix_smoothed(nbins_mean,nbins_mean);
  for (int i=0;i<nbins_mean;i++) {
    for (int j=0;j<nbins_mean;j++) {
      correlation_matrix(i,j)          = (*cov_mat_bootstrapping)(i,j) / sqrt((*cov_mat_bootstrapping)(i,i) * (*cov_mat_bootstrapping)(j,j));
      correlation_matrix_smoothed(i,j) = cov_mat_bootstrapping_temp(i,j)      / sqrt(cov_mat_bootstrapping_temp(i,i)      * cov_mat_bootstrapping_temp(j,j));
      if ((*cov_mat_bootstrapping)(i,j)==0) { correlation_matrix(i,j)     = 0; }
      if (cov_mat_bootstrapping_temp(i,j)==0) { correlation_matrix_smoothed(i,j) = 0; }
      if (i==j) {
        correlation_matrix(i,j) = 1;
        correlation_matrix_smoothed(i,j) = 1;
      }
    }
  }

  std::vector<GPPoint> gp_points_debug;
  double pt_arr_debug[5] = {0,0,0,0,0};
  for (int i=0;i<180;i++) {
    pt_arr_debug[2] = 50 + 10*i;
    gp_points_debug.push_back(GPPoint(pt_arr_debug));
  }

  reg.Predict(gp_points_debug);
  TVectorD vec_mean_debug              = reg.PosteriorMean();
  TMatrixD cov_mat_bootstrapping_debug = reg.PosteriorCov();

  TFile* debug_file = new TFile("/uboone/data/users/lcoopert/LEE/LEEana_xs_1D_Nov12_2021/hist_rootfiles/DetVar/debug_smoothing.root","RECREATE");
  debug_file->cd();
  vec_mean_debug.Write("vec_mean_debug");
  vec_mean->Write("vec_mean_diff");
  vec_mean_temp.Write("vec_mean_diff_smoothed");
  correlation_matrix.Write("correlation_matrix");
  correlation_matrix_smoothed.Write("correlation_matrix_smoothed");
  cov_mat_bootstrapping->Write("cov_mat_bootstrapping");
  cov_mat_bootstrapping_scaled.Write("cov_mat_scaled");
  cov_mat_bootstrapping_debug.Write("cov_mat_smoothed");
  debug_file->Write();
  debug_file->Close();
  //------------
  

  //Store output,  MUST BE LAST
  int nrows = cov_mat_bootstrapping->GetNrows();
  for (int i=0;i<nrows;i++) {
      (*vec_mean)[i] = vec_mean_temp[i];
    for (int j=0;j<nrows;j++) {
      (*cov_mat_bootstrapping)(i,j) = cov_mat_bootstrapping_temp(i,j);
    }
  }
}



