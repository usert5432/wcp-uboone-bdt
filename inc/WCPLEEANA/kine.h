#ifndef UBOONE_LEE_KINE
#define UBOONE_LEE_KINE

#define DLEE_NAME(version, flavor, contain, energy) \
    vlne_ ##version## _ ##flavor## _ ##contain## _ ##energy## E

/* Dirty hack to make string subsitution work */
#define TO_STRING_BASE(x) #x
#define TO_STRING(x) TO_STRING_BASE(x)

#define DLEE_NAME_S(version, flavor, contain, energy) \
    TO_STRING(DLEE_NAME(version, flavor, contain, energy))

#define DLEE_DEFINE(version)                            \
    bool has_vlne_##version;                            \
    Float_t DLEE_NAME(version, numu, full, primary);    \
    Float_t DLEE_NAME(version, numu, full, total);      \
    Float_t DLEE_NAME(version, numu, partial, primary); \
    Float_t DLEE_NAME(version, numu, partial, total);

#define DLEE_INIT(version, info)                                \
    {                                                           \
        info.DLEE_NAME(version, numu, full, primary)    = 0;    \
        info.DLEE_NAME(version, numu, full, total)      = 0;    \
        info.DLEE_NAME(version, numu, partial, primary) = 0;    \
        info.DLEE_NAME(version, numu, partial, total)   = 0;    \
    }

#define DLEE_INFO_PUT_TREE_ADDRESS(version, tree, info)             \
    {                                                               \
        if (info.has_vlne_##version)                                \
        {                                                           \
            tree->Branch(                                           \
                DLEE_NAME_S(version, numu, full, primary),          \
                &info.DLEE_NAME(version, numu, full, primary),      \
                "data/F"                                            \
            );                                                      \
            tree->Branch(                                           \
                DLEE_NAME_S(version, numu, full, total),            \
                &info.DLEE_NAME(version, numu, full, total),        \
                "data/F"                                            \
            );                                                      \
            tree->Branch(                                           \
                DLEE_NAME_S(version, numu, partial, primary),       \
                &info.DLEE_NAME(version, numu, partial, primary),   \
                "data/F"                                            \
            );                                                      \
            tree->Branch(                                           \
                DLEE_NAME_S(version, numu, partial, total),         \
                &info.DLEE_NAME(version, numu, partial, total),     \
                "data/F"                                            \
            );                                                      \
        }                                                           \
    }

#define DLEE_INFO_SET_TREE_ADDRESS(version, tree, info)             \
    {                                                               \
        info.has_vlne_##version = false;                            \
        if (tree->GetBranch(DLEE_NAME_S(version, numu, full, primary))) \
        {                                                           \
            tree->SetBranchAddress(                                 \
                DLEE_NAME_S(version, numu, full, primary),          \
                &info.DLEE_NAME(version, numu, full, primary)       \
            );                                                      \
            tree->SetBranchAddress(                                 \
                DLEE_NAME_S(version, numu, full, total),            \
                &info.DLEE_NAME(version, numu, full, total)         \
            );                                                      \
            tree->SetBranchAddress(                                 \
                DLEE_NAME_S(version, numu, partial, primary),       \
                &info.DLEE_NAME(version, numu, partial, primary)    \
            );                                                      \
            tree->SetBranchAddress(                                 \
                DLEE_NAME_S(version, numu, partial, total),         \
                &info.DLEE_NAME(version, numu, partial, total)      \
            );                                                      \
            info.has_vlne_##version = true;                         \
        }                                                           \
    }

#define DLEE_SET_BRANCH_STATUS(version, tree, status)                   \
    {                                                                   \
        if (tree->GetBranch(DLEE_NAME_S(version, numu, full, primary))) \
        {                                                               \
            tree->SetBranchStatus(                                      \
                DLEE_NAME_S(version, numu, full, primary), status       \
            );                                                          \
            tree->SetBranchStatus(                                      \
                DLEE_NAME_S(version, numu, full, total), status         \
            );                                                          \
            tree->SetBranchStatus(                                      \
                DLEE_NAME_S(version, numu, partial, primary), status    \
            );                                                          \
            tree->SetBranchStatus(                                      \
                DLEE_NAME_S(version, numu, partial, total), status      \
            );                                                          \
        }                                                               \
    }


namespace LEEana{
struct KineInfo{
  Float_t kine_reco_Enu;
  Float_t kine_reco_add_energy;
  std::vector<float> *kine_energy_particle;
  std::vector<int> *kine_energy_info; 
  std::vector<int> *kine_particle_type;
  std::vector<int> *kine_energy_included;
  Float_t kine_pio_mass;
  Int_t kine_pio_flag;
  Float_t kine_pio_vtx_dis;
  Float_t kine_pio_energy_1;
  Float_t kine_pio_theta_1;
  Float_t kine_pio_phi_1;
  Float_t kine_pio_dis_1;
  Float_t kine_pio_energy_2;
  Float_t kine_pio_theta_2;
  Float_t kine_pio_phi_2;
  Float_t kine_pio_dis_2;
  Float_t kine_pio_angle;

  // DL energy
  bool has_dl_ee;
  DLEE_DEFINE(v1)
  DLEE_DEFINE(v2)
  DLEE_DEFINE(v3)
  DLEE_DEFINE(v4)
  DLEE_DEFINE(v5)
  DLEE_DEFINE(v6)
  DLEE_DEFINE(v7)
};

void set_tree_address(TTree *tree0, KineInfo& tagger_info);
void put_tree_address(TTree *Tsig, KineInfo& tagger_info);
 void clear_kine_info(KineInfo& tagger_info);
}

void LEEana::clear_kine_info(KineInfo& tagger_info){
  tagger_info.kine_reco_Enu=0;
  tagger_info.kine_reco_add_energy=0;
  tagger_info.kine_energy_particle->clear();
  tagger_info.kine_energy_info->clear(); 
  tagger_info.kine_particle_type->clear();
  tagger_info.kine_energy_included->clear();
  tagger_info.kine_pio_mass=0;
  tagger_info.kine_pio_flag=0;
  tagger_info.kine_pio_vtx_dis=0;
  tagger_info.kine_pio_energy_1=0;
  tagger_info.kine_pio_theta_1=0;
  tagger_info.kine_pio_phi_1=0;
  tagger_info.kine_pio_dis_1=0;
  tagger_info.kine_pio_energy_2=0;
  tagger_info.kine_pio_theta_2=0;
  tagger_info.kine_pio_phi_2=0;
  tagger_info.kine_pio_dis_2=0;
  tagger_info.kine_pio_angle=0;

  DLEE_INIT(v1, tagger_info);
  DLEE_INIT(v2, tagger_info);
  DLEE_INIT(v3, tagger_info);
  DLEE_INIT(v4, tagger_info);
  DLEE_INIT(v5, tagger_info);
  DLEE_INIT(v6, tagger_info);
  DLEE_INIT(v7, tagger_info);
}

void LEEana::set_tree_address(TTree *tree0, KineInfo& tagger_info){
  tree0->SetBranchAddress("kine_reco_Enu", &tagger_info.kine_reco_Enu);
  tree0->SetBranchAddress("kine_reco_add_energy", &tagger_info.kine_reco_add_energy);
  tree0->SetBranchAddress("kine_energy_particle", &tagger_info.kine_energy_particle);
  tree0->SetBranchAddress("kine_energy_info", &tagger_info.kine_energy_info); 
  tree0->SetBranchAddress("kine_particle_type", &tagger_info.kine_particle_type);
  tree0->SetBranchAddress("kine_energy_included", &tagger_info.kine_energy_included);
  tree0->SetBranchAddress("kine_pio_mass", &tagger_info.kine_pio_mass);
  tree0->SetBranchAddress("kine_pio_flag", &tagger_info.kine_pio_flag);
  tree0->SetBranchAddress("kine_pio_vtx_dis", &tagger_info.kine_pio_vtx_dis);
  tree0->SetBranchAddress("kine_pio_energy_1", &tagger_info.kine_pio_energy_1);
  tree0->SetBranchAddress("kine_pio_theta_1", &tagger_info.kine_pio_theta_1);
  tree0->SetBranchAddress("kine_pio_phi_1", &tagger_info.kine_pio_phi_1);
  tree0->SetBranchAddress("kine_pio_dis_1", &tagger_info.kine_pio_dis_1);
  tree0->SetBranchAddress("kine_pio_energy_2", &tagger_info.kine_pio_energy_2);
  tree0->SetBranchAddress("kine_pio_theta_2", &tagger_info.kine_pio_theta_2);
  tree0->SetBranchAddress("kine_pio_phi_2", &tagger_info.kine_pio_phi_2);
  tree0->SetBranchAddress("kine_pio_dis_2", &tagger_info.kine_pio_dis_2);
  tree0->SetBranchAddress("kine_pio_angle", &tagger_info.kine_pio_angle);

    DLEE_INFO_SET_TREE_ADDRESS(v1, tree0, tagger_info);
    DLEE_INFO_SET_TREE_ADDRESS(v2, tree0, tagger_info);
    DLEE_INFO_SET_TREE_ADDRESS(v3, tree0, tagger_info);
    DLEE_INFO_SET_TREE_ADDRESS(v4, tree0, tagger_info);
    DLEE_INFO_SET_TREE_ADDRESS(v5, tree0, tagger_info);
    DLEE_INFO_SET_TREE_ADDRESS(v6, tree0, tagger_info);
    DLEE_INFO_SET_TREE_ADDRESS(v7, tree0, tagger_info);
}

void LEEana::put_tree_address(TTree *tree0, KineInfo& tagger_info){
  tree0->Branch("kine_reco_Enu", &tagger_info.kine_reco_Enu,"data/F");
  tree0->Branch("kine_reco_add_energy", &tagger_info.kine_reco_add_energy,"data/F");
  tree0->Branch("kine_energy_particle", &tagger_info.kine_energy_particle);
  tree0->Branch("kine_energy_info", &tagger_info.kine_energy_info); 
  tree0->Branch("kine_particle_type", &tagger_info.kine_particle_type);
  tree0->Branch("kine_energy_included", &tagger_info.kine_energy_included);
  tree0->Branch("kine_pio_mass", &tagger_info.kine_pio_mass,"data/F");
  tree0->Branch("kine_pio_flag", &tagger_info.kine_pio_flag,"data/I");
  tree0->Branch("kine_pio_vtx_dis", &tagger_info.kine_pio_vtx_dis,"data/F");
  tree0->Branch("kine_pio_energy_1", &tagger_info.kine_pio_energy_1,"data/F");
  tree0->Branch("kine_pio_theta_1", &tagger_info.kine_pio_theta_1,"data/F");
  tree0->Branch("kine_pio_phi_1", &tagger_info.kine_pio_phi_1,"data/F");
  tree0->Branch("kine_pio_dis_1", &tagger_info.kine_pio_dis_1,"data/F");
  tree0->Branch("kine_pio_energy_2", &tagger_info.kine_pio_energy_2,"data/F");
  tree0->Branch("kine_pio_theta_2", &tagger_info.kine_pio_theta_2,"data/F");
  tree0->Branch("kine_pio_phi_2", &tagger_info.kine_pio_phi_2,"data/F");
  tree0->Branch("kine_pio_dis_2", &tagger_info.kine_pio_dis_2,"data/F");
  tree0->Branch("kine_pio_angle", &tagger_info.kine_pio_angle,"data/F");

    DLEE_INFO_PUT_TREE_ADDRESS(v1, tree0, tagger_info);
    DLEE_INFO_PUT_TREE_ADDRESS(v2, tree0, tagger_info);
    DLEE_INFO_PUT_TREE_ADDRESS(v3, tree0, tagger_info);
    DLEE_INFO_PUT_TREE_ADDRESS(v4, tree0, tagger_info);
    DLEE_INFO_PUT_TREE_ADDRESS(v5, tree0, tagger_info);
    DLEE_INFO_PUT_TREE_ADDRESS(v6, tree0, tagger_info);
    DLEE_INFO_PUT_TREE_ADDRESS(v7, tree0, tagger_info);
}


#endif
