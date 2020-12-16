#ifndef UBOONE_LEE_CUTS
#define UBOONE_LEE_CUTS

// define cuts here ...
#include "TCut.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "tagger.h"
#include "kine.h"
#include "eval.h"
#include "pfeval.h"

#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

namespace LEEana{

  double get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, TString var_name="kine_reco_Enu");
  bool get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  double get_weight(TString weight_name, EvalInfo& eval);

  int get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  
  // generic neutrino cuts
  // TCut generic_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length >15";
  bool is_generic(EvalInfo& info);
  
  // preselection cuts
  // TCut preselect_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length > 0";
  bool is_preselection(EvalInfo& info); 
  
  // nueCC cuts
  // TCut nueCC_cut = "numu_cc_flag >=0 && nue_score > 7.0";
  bool is_nueCC(TaggerInfo& tagger_info);
  bool is_loosenueCC(TaggerInfo& tagger_info);
  
  bool is_far_sideband(KineInfo& kine, TaggerInfo& tagger);
  bool is_near_sideband(KineInfo& kine, TaggerInfo& tagger);
  
  // numuCC cuts
  // TCut numuCC_cut = "numu_cc_flag >=0 && numu_score > 0.9";
  bool is_numuCC(TaggerInfo& tagger_info);
  bool is_numuCC_tight(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_numuCC_1mu0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_numuCC_lowEhad(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_numuCC_cutbased(TaggerInfo& tagger_info);
  
  // pio cuts (with and without vertex)
  // TCut pi0_cut = "(kine_pio_flag==1 && kine_pio_vtx_dis < 9 || kine_pio_flag ==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0  && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
  bool is_pi0(KineInfo& kine);
  
  // must be with vertex ...
  // TCut cc_pi0_cut = "(kine_pio_flag==1 && kine_pio_vtx_dis < 9 || kine_pio_flag ==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0  && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
  bool is_cc_pi0(KineInfo& kine);
  
  

 
  // NC cuts
  // TCut NC_cut = "(!cosmict_flag) && numu_score < 0.0";
  bool is_NC(TaggerInfo& tagger_info);
  
  
  // TCut FC_cut = "match_isFC==1";
  // TCut PC_cut = "match_isFC==0";
  
  bool is_FC(EvalInfo& eval);
  
  
  // TCut truth_nueCC_inside = "abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1";
  // TCut truth_numuCC_inside = "abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1";
  bool is_truth_nueCC_inside(EvalInfo& eval);
  bool is_truth_numuCC_inside(EvalInfo& eval);
}

double LEEana::get_weight(TString weight_name, EvalInfo& eval){
  if (weight_name == "cv_spline"){
    return eval.weight_cv * eval.weight_spline;
  }else if (weight_name == "cv_spline_cv_spline"){
    return pow(eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "unity" || weight_name == "unity_unity"){
    return 1;
  }else if (weight_name == "lee_cv_spline"){
    return (eval.weight_lee * eval.weight_cv * eval.weight_spline);
  }else if (weight_name == "lee_cv_spline_lee_cv_spline"){
    return pow(eval.weight_lee * eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "lee_cv_spline_cv_spline" || weight_name == "cv_spline_lee_cv_spline"){
    return eval.weight_lee * pow(eval.weight_cv * eval.weight_spline,2);
  }else{
    std::cout <<"Unknown weights: " << weight_name << std::endl;
  }
	    
  
  return 1;
}

double LEEana::get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, TString var_name ){
  if (var_name == "kine_reco_Enu"){
    return kine.kine_reco_Enu;
  }else if (var_name == "match_energy"){
    return eval.match_energy;
  }else if (var_name == "pi0_energy"){
    double pi0_mass = 135;
    double alpha = fabs(kine.kine_pio_energy_1 - kine.kine_pio_energy_2)/(kine.kine_pio_energy_1 + kine.kine_pio_energy_2);
    return pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(kine.kine_pio_angle/180.*3.1415926)))-1);
  }else if (var_name == "pi0_mass"){
    return kine.kine_pio_mass;
  }else if (var_name == "nue_score"){
    return tagger.nue_score;
  }else if (var_name == "numu_score"){
    return tagger.numu_score;
  }else if (var_name == "shower_energy"){
    return tagger.mip_energy;
  }else if (var_name == "shower_angle_beam"){
    return tagger.mip_angle_beam;
  }else if (var_name == "shower_angle_vertical"){
    return tagger.spt_angle_vertical;
  }else if (var_name == "shwvtx_nuvtx_dis"){
    return sqrt(pow(pfeval.reco_nuvtxX-pfeval.reco_showervtxX,2)+pow(pfeval.reco_nuvtxY-pfeval.reco_showervtxY,2)+pow(pfeval.reco_nuvtxZ-pfeval.reco_showervtxZ,2)); 
  }else if (var_name == "median_dQdx"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.mip_vec_dQ_dx_2);
    dqdx.push_back(tagger.mip_vec_dQ_dx_3);
    dqdx.push_back(tagger.mip_vec_dQ_dx_4);
    dqdx.push_back(tagger.mip_vec_dQ_dx_5);
    dqdx.push_back(tagger.mip_vec_dQ_dx_6);
    dqdx.push_back(tagger.mip_vec_dQ_dx_7);
    dqdx.push_back(tagger.mip_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    return vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
  }else if (var_name == "median_dEdx"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.mip_vec_dQ_dx_2);
    dqdx.push_back(tagger.mip_vec_dQ_dx_3);
    dqdx.push_back(tagger.mip_vec_dQ_dx_4);
    dqdx.push_back(tagger.mip_vec_dQ_dx_5);
    dqdx.push_back(tagger.mip_vec_dQ_dx_6);
    dqdx.push_back(tagger.mip_vec_dQ_dx_7);
    dqdx.push_back(tagger.mip_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    float median_dqdx = vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
    float alpha = 1.;
    float beta = 0.255;
    float median_dedx = (exp((median_dqdx*43e3) * 23.6e-6*beta/1.38/0.273) - alpha)/(beta/1.38/0.273);
    if(median_dedx<0) median_dedx = 0;
    if(median_dedx>50) median_dedx = 50;
    return median_dedx; // MeV/cm
  }else if (var_name == "reco_showervtxX"){
      return pfeval.reco_showervtxX;
  }else if (var_name == "reco_nuvtxX"){
      return pfeval.reco_nuvtxX;
  }else if (var_name == "reco_nuvtxY"){
      return pfeval.reco_nuvtxY;
  }else if (var_name == "reco_nuvtxZ"){
      return pfeval.reco_nuvtxZ;
  }else if (var_name == "mip_quality_n_tracks"){
      return tagger.mip_quality_n_tracks;
  }else if (var_name == "mip_quality_n_showers"){
      return tagger.mip_quality_n_showers;
  }else if (var_name == "gap_n_bad"){
      return tagger.gap_n_bad;
  }else if (var_name == "muon_KE"){
      return pfeval.reco_muonMomentum[3]*1000.-105.66; // GeV --> MeV
  }else if (var_name == "muon_costheta"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      return TMath::Cos(muonMomentum.Theta());
  }else if (var_name == "muon_theta"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      return muonMomentum.Theta()*180./TMath::Pi();
  }else if (var_name == "muon_phi"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      return muonMomentum.Phi()/TMath::Pi()*180.;
  }else if (var_name == "proton_KE"){
      return pfeval.reco_protonMomentum[3]*1000.-938.27; // GeV--> MeV
  }else if (var_name == "proton_theta"){
      TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
      return protonMomentum.Theta()/TMath::Pi()*180.;
  }else if (var_name == "proton_phi"){
      TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
      return protonMomentum.Phi()/TMath::Pi()*180.;
  }else if (var_name == "Ehadron"){
      /* Float_t Ehadron = kine.kine_reco_Enu; */
      /* for(size_t i=0; i<kine.kine_energy_particle->size(); i++) */
      /* { */
      /*     int pdgcode = kine.kine_particle_type->at(i); */
      /*     if(abs(pdgcode)==13) Ehadron = Ehadron - kine.kine_energy_particle->at(i) - 105.658; */ 
      /*     //if(abs(pdgcode)==11) Ehadron = Ehadron - kine.kine_energy_particle->at(i); */ 
      /* } */
      return kine.kine_reco_Enu - pfeval.reco_muonMomentum[3]*1000.;
  }else if (var_name == "Q2"){
      Float_t Enu = kine.kine_reco_Enu;
      Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
      Float_t Ehadron = Enu - Emu;
      Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
      return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(1000.*1000.); // GeV^2
  }else if (var_name == "x_Bjorken"){
      Float_t Enu = kine.kine_reco_Enu;
      Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
      Float_t Ehadron = Enu - Emu;
      Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
      return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(2*938.272*Ehadron);
  }else if (var_name == "N_tracks"){
      int N_tracks = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==11) continue;
          if(kine.kine_energy_particle->at(i)<10) continue;
          if(abs(pdgcode)==13 || abs(pdgcode)==211){
            N_tracks += 1;
          }
          else if(kine.kine_energy_particle->at(i)>35){ // proton KE threshold
              N_tracks += 1; 
          }
      }
      return N_tracks;
  }else if (var_name == "N_showers"){
      int N_showers = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)!=11) continue;
          if(kine.kine_energy_particle->at(i)>10) N_showers += 1;
      }
      return N_showers;
  }else{
    std::cout << "No such variable: " << var_name << std::endl;
    exit(EXIT_FAILURE);
  }
  return -1;
}

int LEEana::get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){
  for (auto it = map_cut_xs_bin.begin(); it != map_cut_xs_bin.end(); it++){
    TString cut_name = it->first;
    int number = it->second;

    if (cut_file == 1){
      if (cut_name == "numuCC.inside.Enu.le.300"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=300) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.400"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=400 ) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.400.gt.300"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=400 && eval.truth_nuEnergy>300) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.500.gt.400"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=500 && eval.truth_nuEnergy>400) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.600.gt.500"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=600 && eval.truth_nuEnergy>500) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.700.gt.600"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=700 && eval.truth_nuEnergy>600) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.800.gt.700"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=800 && eval.truth_nuEnergy>700) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.900.gt.800"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=900 && eval.truth_nuEnergy>800) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1000.gt.900"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1000 && eval.truth_nuEnergy>900) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1100.gt.1000"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1100 && eval.truth_nuEnergy>1000) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1100"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1100) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1000"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1000) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1500.gt.1200"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1500 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1400.gt.1200"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1400 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1600.gt.1400"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1600 && eval.truth_nuEnergy>1400) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2000.gt.1600"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2000 && eval.truth_nuEnergy>1600) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2500.gt.2000"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2500 && eval.truth_nuEnergy>2000) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.2500"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2500) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.1500"){
	if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>1500) return number;
      }else{
	std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    
  }
  
  return -1;
}

bool LEEana::get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){

  bool flag_truth_inside = false; // in the active volume
  if (eval.truth_vtxX > -1 && eval.truth_vtxX <= 254.3 &&  eval.truth_vtxY >-115.0 && eval.truth_vtxY<=117.0 && eval.truth_vtxZ > 0.6 && eval.truth_vtxZ <=1036.4) flag_truth_inside = true;

  // definition of additional cuts
  std::map<std::string, bool> map_cuts_flag;
  if(is_far_sideband(kine, tagger)) map_cuts_flag["farsideband"] = true; 
  else map_cuts_flag["farsideband"] = false; 
  
  if(is_near_sideband(kine, tagger)) map_cuts_flag["nearsideband"] = true; 
  else map_cuts_flag["nearsideband"] = false; 
 
  if(is_nueCC(tagger)) map_cuts_flag["nueCC"] = true;
  else map_cuts_flag["nueCC"] = false;

  if(is_loosenueCC(tagger)) map_cuts_flag["loosenueCC"] = true;
  else map_cuts_flag["loosenueCC"] = false;
  
  if(is_generic(eval)) map_cuts_flag["generic"] = true;
  else map_cuts_flag["generic"] = false;

  if(eval.truth_nuEnergy <=400) map_cuts_flag["LowEintnueCC"] = true;
  else map_cuts_flag["LowEintnueCC"] = false;
  
  if (!(eval.truth_nuEnergy <=400)) map_cuts_flag["antiLowEintnueCC"] = true;
  else map_cuts_flag["antiLowEintnueCC"] = false;

  if(eval.truth_nuEnergy<=400) map_cuts_flag["LowEnu"] = true;
  else map_cuts_flag["LowEnu"] = false;
  
  if(!(eval.truth_nuEnergy<=400)) map_cuts_flag["antiLowEnu"] = true;
  else map_cuts_flag["antiLowEnu"] = false;

  if(eval.match_completeness_energy/eval.truth_energyInside<0.1) map_cuts_flag["badmatch"] = true;
  else map_cuts_flag["badmatch"] = false;
  
  if(eval.match_completeness_energy/eval.truth_energyInside>=0.1 && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["numuCCinFV"] = true;
  else map_cuts_flag["numuCCinFV"] = false;

   if(eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 ) map_cuts_flag["XsnumuCCinFV"] = true;
  else map_cuts_flag["XsnumuCCinFV"] = false;

  if(eval.match_completeness_energy/eval.truth_energyInside>=0.1 && abs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["nueCCinFV"] = true;
  else map_cuts_flag["nueCCinFV"] = false;
    
  if(eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["NCinFV"] = true;
  else map_cuts_flag["NCinFV"] = false;

  if(eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_vtxInside==0) map_cuts_flag["outFV"] = true;
  else map_cuts_flag["outFV"] = false;
      
  if(eval.match_completeness_energy/eval.truth_energyInside>=0.1 && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0) map_cuts_flag["CCpi0inFV"] = true;
  else map_cuts_flag["CCpi0inFV"] = false;
      
  if (eval.match_completeness_energy/eval.truth_energyInside>=0.1 && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0) map_cuts_flag["NCpi0inFV"] = true;
  else map_cuts_flag["NCpi0inFV"] = false;

  if(pfeval.truth_nuIntType == 1000) map_cuts_flag["MEC"] = true; // 10% of this type from unclassified res?
  else map_cuts_flag["MEC"] = false;
  
  if(pfeval.truth_nuIntType == 1001) map_cuts_flag["CCQE"] = true;
  else map_cuts_flag["CCQE"] = false;

  if(pfeval.truth_nuIntType == 1002) map_cuts_flag["NCQE"] = true;
  else map_cuts_flag["NCQE"] = false;

  if(pfeval.truth_nuIntType >= 1003 && pfeval.truth_nuIntType<=1090 && !(pfeval.truth_nuIntType>=1006 && pfeval.truth_nuIntType<=1009) && !(pfeval.truth_nuIntType>=1013 && pfeval.truth_nuIntType<=1016)) map_cuts_flag["CCRES"] = true;
  else map_cuts_flag["CCRES"] = false;

  if((pfeval.truth_nuIntType>=1006 && pfeval.truth_nuIntType<=1009) || (pfeval.truth_nuIntType>=1013 && pfeval.truth_nuIntType<=1016)) map_cuts_flag["NCRES"] = true;
  else map_cuts_flag["NCRES"] = false;

  if(pfeval.truth_nuIntType == 1091 || pfeval.truth_nuIntType == 1092) map_cuts_flag["DIS"] = true;
  else map_cuts_flag["DIS"] = false;

  if(pfeval.truth_nuIntType<1000 || pfeval.truth_nuIntType>1092) map_cuts_flag["otherXs"] = true;
  else map_cuts_flag["otherXs"] = false;


  // figure out additional cuts and flag_data ...
  bool flag_add = true;
  if(add_cut == "all") flag_add = true;
  else if( (flag_data && (add_cut=="farsideband" || add_cut=="nearsideband" || add_cut=="nueCC" || add_cut=="generic" || add_cut=="loosenueCC")) || !flag_data ){ 
      std::istringstream sss(add_cut.Data());
      for(std::string line; std::getline(sss, line, '_');){
          if(map_cuts_flag.find(line)!=map_cuts_flag.end()){
              flag_add *= map_cuts_flag[line];
          }
          else{
              std::cout<<"ERROR: add_cut "<<line<<" not defined!\n";
              exit(EXIT_FAILURE);
          }
      } 
  }
  else{ 
    std::cout<<"ERROR: add_cut "<<add_cut<<" of channel "<< ch_name <<" is not assigned to sample "<<flag_data<<" [1: data; 0: mc]\n";
    std::cout<<"Please modify inc/WCPLEEANA/cuts.h\n";
    exit(EXIT_FAILURE);
  }
 
  if (!flag_add) return false;

  bool flag_generic = is_generic(eval);
  bool flag_numuCC = is_numuCC(tagger);
  bool flag_numuCC_tight = is_numuCC_tight(tagger, pfeval);
  bool flag_numuCC_1mu0p = is_numuCC_1mu0p(tagger, kine, pfeval);
  bool flag_numuCC_lowEhad = is_numuCC_lowEhad(tagger, kine, pfeval);
  bool flag_numuCC_cutbased = is_numuCC_cutbased(tagger);
  bool flag_nueCC = is_nueCC(tagger);
  bool flag_nueCC_loose = is_loosenueCC(tagger);
  bool flag_pi0 = is_pi0(kine);
  bool flag_cc_pi0 = is_cc_pi0(kine);
  bool flag_NC = is_NC(tagger);
  bool flag_FC = is_FC(eval);

  if (ch_name == "LEE_FC_nueoverlay"  || ch_name == "nueCC_FC_nueoverlay"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_ext" || ch_name == "BG_nueCC_FC_dirt" || ch_name =="nueCC_FC_bnb"){
    //nueCC FC
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_overlay"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "LEE_PC_nueoverlay" || ch_name == "nueCC_PC_nueoverlay" ){
    // nueCC PC
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_ext" || ch_name == "BG_nueCC_PC_dirt" || ch_name == "nueCC_PC_bnb"){
    // nueCC PC
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_overlay"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_FC_overlay" || ch_name =="BG_CCpi0_nonueCC_FC_ext" || ch_name == "BG_CCpi0_nonueCC_FC_dirt" || ch_name == "CCpi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_PC_overlay" || ch_name == "BG_CCpi0_nonueCC_PC_ext" || ch_name == "BG_CCpi0_nonueCC_PC_dirt" || ch_name == "CCpi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "NCpi0_nonueCC_overlay" || ch_name == "BG_NCpi0_nonueCC_ext" || ch_name == "BG_NCpi0_nonueCC_dirt" || ch_name == "NCpi0_nonueCC_bnb"){
    if (flag_NC && flag_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "nueCC_bnb" || ch_name == "nueCC_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "nueCC_nueoverlay" || ch_name == "nueCC_bnb") return true;
    else return false;
  }else if (ch_name == "all_but_nueCC_bnb" || ch_name == "all_but_nueCC_overlay" || ch_name == "all_but_nueCC_ext" || ch_name == "all_but_nueCC_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "all_but_nueCC_overlay" || ch_name != "all_but_nueCC_overlay") return true;
    else return false;
  }else if (ch_name == "nueCC_bnb1" || ch_name == "nueCC_nueoverlay1"){
    if (flag_truth_inside &&  ch_name == "nueCC_nueoverlay1" || ch_name == "nueCC_bnb1") return true;
    else return false;
  }else if (ch_name == "all_but_nueCC_bnb1" || ch_name == "all_but_nueCC_overlay1" || ch_name == "all_but_nueCC_ext1" || ch_name == "all_but_nueCC_dirt1"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "all_but_nueCC_overlay1" || ch_name != "all_but_nueCC_overlay1") return true;
    else return false;
  }else if (ch_name == "testA_bnb" || ch_name == "testA_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testA_nueoverlay" || ch_name == "testA_bnb") return true;
    else return false;
  }else if (ch_name == "testA_overlay" || ch_name == "testA_ext" || ch_name == "testA_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testA_overlay" || ch_name != "testA_overlay") return true;
    else return false;
  }else if (ch_name == "testB_bnb" || ch_name == "testB_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testB_nueoverlay" || ch_name == "testB_bnb") return true;
    else return false;
  }else if (ch_name == "testB_overlay" || ch_name == "testB_ext" || ch_name == "testB_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testB_overlay" || ch_name != "testB_overlay") return true;
    else return false;
  }else if (ch_name == "testC_bnb" || ch_name == "testC_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testC_nueoverlay" || ch_name == "testC_bnb") return true;
    else return false;
  }else if (ch_name == "testC_overlay" || ch_name == "testC_ext" || ch_name == "testC_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testC_overlay" || ch_name != "testC_overlay") return true;
    else return false;
  }else if (ch_name == "testD_bnb" || ch_name == "testD_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testD_nueoverlay" || ch_name == "testD_bnb") return true;
    else return false;
  }else if (ch_name == "testD_overlay" || ch_name == "testD_ext" || ch_name == "testD_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testD_overlay" || ch_name != "testD_overlay") return true;
    else return false;
 // Janet's requests: <600 MeV numuCC PC, FC for three variables = 6 obs channels 
  }else if (ch_name == "numuCC_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=0 && kine.kine_reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=0 && kine.kine_reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC2_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=0 && kine.kine_reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC2_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=0 && kine.kine_reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC3_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=0 && kine.kine_reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC3_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=0 && kine.kine_reco_Enu<600) return true;
    else return false;

  }else if (ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=600 && kine.kine_reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=600 && kine.kine_reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=600 && kine.kine_reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=600 && kine.kine_reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=600 && kine.kine_reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && kine.kine_reco_Enu>=600 && kine.kine_reco_Enu<1500) return true;
    else return false;

  }else if (ch_name == "numuCC_extra_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_extra_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_extra_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_extra_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_extra_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_extra_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_extra_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC_extra2_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_extra2_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_extra2_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_extra2_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra2_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_extra2_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_extra2_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_extra2_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC2_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC2_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC3_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC3_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC4_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC4_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC4_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC4_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC4_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC4_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC4_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC4_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC2_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC2_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC3_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC3_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC4_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC4_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC4_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC4_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC4_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC4_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC4_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC4_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

 // Mike Shaevitz >800 MeV nueCC PC+FC 1 obs channel
  }else if (ch_name == "nueCC_extra_nueoverlay"){
    if (flag_nueCC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_ext" || ch_name == "BG_nueCC_extra_dirt" || ch_name =="nueCC_extra_bnb"){
    if (flag_nueCC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_overlay"){
    if (flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
 // cut-based numuCC FC/PC 2 obs channels   
  }else if (ch_name == "numuCC_cutbased_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_cutbased_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_cutbased_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_cutbased_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_cutbased && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_cutbased_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_cutbased_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_cutbased_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_cutbased_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_cutbased && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

 // generic selection nu PC+FC 1 obs channel   
  }else if (ch_name == "generic_nu_overlay" || ch_name == "BG_generic_nu_ext" || ch_name =="BG_generic_nu_dirt" || ch_name == "generic_nu_bnb"){
    if (flag_generic) return true;
    else return false;
 // numuCC selection PC+FC 1 obs channel   
  }else if (ch_name == "numuCC_overlay" || ch_name == "BG_numuCC_ext" || ch_name =="BG_numuCC_dirt" || ch_name == "numuCC_bnb"){
    if (flag_numuCC) return true;
    else return false;
 // cutbased numuCC selection PC+FC 1 obs channel   
  }else if (ch_name == "numuCC_cutbased_overlay" || ch_name == "BG_numuCC_cutbased_ext" || ch_name =="BG_numuCC_cutbased_dirt" || ch_name == "numuCC_cutbased_bnb"){
    if (flag_numuCC_cutbased) return true;
    else return false;
 // nueCC 3 variables: n_trakcs, n_showers, gap_n_bad, FC/PC x3 = 6 channels; 4 additional channels 
  }else if (ch_name == "nueCC2_FC_nueoverlay"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_FC_ext" || ch_name == "BG_nueCC2_FC_dirt" || ch_name =="nueCC2_FC_bnb"){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_FC_overlay"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "nueCC2_PC_nueoverlay" ){
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_PC_ext" || ch_name == "BG_nueCC2_PC_dirt" || ch_name == "nueCC2_PC_bnb"){
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_PC_overlay"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "nueCC3_FC_nueoverlay"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_FC_ext" || ch_name == "BG_nueCC3_FC_dirt" || ch_name =="nueCC3_FC_bnb"){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_FC_overlay"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "nueCC3_PC_nueoverlay" ){
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_PC_ext" || ch_name == "BG_nueCC3_PC_dirt" || ch_name == "nueCC3_PC_bnb"){
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_PC_overlay"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
    // add some cuts for Xs related cases ...
  }else if (ch_name == "numuCC_FC_bnb" || ch_name == "BG_numuCC_FC_ext" || ch_name == "BG_numuCC_FC_dirt"){
    if (flag_numuCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "numuCC_PC_bnb" || ch_name == "BG_numuCC_PC_ext" || ch_name == "BG_numuCC_PC_dirt"){
    if (flag_numuCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "numuCC_signal_FC_overlay" || ch_name == "numuCC_signal_PC_overlay" || ch_name == "numuCC_background_FC_overlay" || ch_name == "numuCC_background_PC_overlay"){
    if (ch_name == "numuCC_signal_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["XsnumuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["XsnumuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["XsnumuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["XsnumuCCinFV"])) return true;
    }
    return false;
  }else{
    std::cout << "Not sure what cut: " << ch_name << std::endl;
  }
  
  return false;
}


bool LEEana::is_far_sideband(KineInfo& kine, TaggerInfo& tagger){
  bool flag = false;

  bool flag_numuCC = is_numuCC(tagger);
  bool flag_pi0 = is_pi0(kine);
  bool flag_cc_pi0 = is_cc_pi0(kine);
  bool flag_NC = is_NC(tagger);

  if ((kine.kine_reco_Enu>=800 && tagger.nue_score >=0) ||
      (tagger.nue_score<=0 && (flag_numuCC || (flag_pi0 && flag_NC) ))) flag = true;
  return flag;
}
bool LEEana::is_near_sideband(KineInfo& kine, TaggerInfo& tagger){
  bool flag = false;
  if (kine.kine_reco_Enu < 800 && tagger.nue_score>0 && (kine.kine_reco_Enu>=600 || tagger.nue_score<=7)) flag = true;
  
  return flag ;
}




bool LEEana::is_truth_nueCC_inside(EvalInfo& eval){
  bool flag = false;

  if (fabs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1)
    flag = true;
  
  return flag;
}

bool LEEana::is_truth_numuCC_inside(EvalInfo& eval){
   bool flag = false;

  if (fabs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1)
    flag = true;
  
  return flag;
}



bool LEEana::is_FC(EvalInfo& eval){
  if (eval.match_isFC){
    return true;
  }else{
    return false;
  }
}

bool LEEana::is_cc_pi0(KineInfo& kine){
  bool flag = false;

  if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300)
    flag = true;

  
  return flag;
}


bool LEEana::is_pi0(KineInfo& kine){
  bool flag = false;

  if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 || kine.kine_pio_flag==2) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300)
    flag = true;

  
  return flag;
}



bool LEEana::is_NC(TaggerInfo& tagger_info){
  bool flag = false;
  if ((!tagger_info.cosmict_flag) && tagger_info.numu_score < 0)
    flag = true;
  
  return flag;
}


bool LEEana::is_numuCC(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9)
    flag = true;
  
  return flag;
}

bool LEEana::is_numuCC_tight(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0)
    flag = true;
  
  return flag;
}

bool LEEana::is_numuCC_1mu0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;
  
  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0){ 
      // 1 lepton <=1 proton 0 charged pion
      // 1 lepton guaranteed by numu cc flag
      // using pi0 flag to remove pi0 component in channel definition
      int Nproton = 0;
      int Npion = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm? 
          if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV 
      }
      if(Nproton==0) flag = true;
  } 
  
  return flag;
}


bool LEEana::is_numuCC_lowEhad(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
    bool flag = false;
   
    if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0){ 
        Float_t Ehadron = kine.kine_reco_Enu - pfeval.reco_muonMomentum[3]*1000.;
        if(Ehadron<200) // MeV
        {
            flag = true;
        }
    }
    return flag;
}

bool LEEana::is_numuCC_cutbased(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag==1 && tagger_info.cosmict_flag==0) 
    flag = true;
  
  return flag;
}


bool LEEana::is_nueCC(TaggerInfo& tagger_info){
  bool flag = false;
  // default 7.0
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 7.0)
    //  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score <= 7.0 && tagger_info.nue_score > 0)
    flag = true;
  
  return flag;
}

bool LEEana::is_loosenueCC(TaggerInfo& tagger_info){
  bool flag = false;
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 4.0)
    flag = true;
  
  return flag;
}

bool LEEana::is_generic(EvalInfo& eval){
  // not very useful for the main analysis
  bool flag = is_preselection(eval);

  flag = flag && (eval.stm_clusterlength > 15);
  return flag;
}

bool LEEana::is_preselection(EvalInfo& eval){
  bool flag = false;

  // match code ...
  int tmp_match_found = eval.match_found;
  if (eval.is_match_found_int){
    tmp_match_found = eval.match_found_asInt;
  }

  if (tmp_match_found == 1 && eval.stm_eventtype != 0 && eval.stm_lowenergy ==0 && eval.stm_LM ==0 && eval.stm_TGM ==0 && eval.stm_STM==0 && eval.stm_FullDead == 0 && eval.stm_clusterlength >0) flag = true;
  
  
  return flag;
}


#endif
