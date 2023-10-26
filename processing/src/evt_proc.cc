#include "cms_hh_proc_interface/processing/interface/evt_proc.hh"

EvtProc::EvtProc(bool return_all, std::vector<std::string> requested, bool use_deep_bjet_wps) {
    _all = return_all;
    _requested = requested;
    _feat_comp = new FeatComp(return_all, _requested, use_deep_bjet_wps);
}

EvtProc::~EvtProc() {}

std::map<std::string, float> EvtProc::process(const LorentzVector& b_1,
                                              const LorentzVector& b_2,
                                              const LorentzVector& l_1,
                                              const LorentzVector& l_2,
                                              const LorentzVector& met,
                                              const LorentzVector& svfit,
                                              const LorentzVector& vbf_1,
                                              const LorentzVector& vbf_2,
                                              const LorentzVector& Nu_1,
                                              const LorentzVector& Nu_2,
                                              const float& HHKin_mass_rawass,
                                              const float& hh_kinfit_chi2,
                                              const float& mt2,
                                              const bool&  is_boosted,
                                              Year year,
                                              const float& res_mass,
                                              Spin spin,
                                              const int& n_vbf,
                                              const int& pairType,
                                              const int& dau1_decayMode,
                                              const int& dau2_decayMode,
                                              const int& dau1_flav,
                                              const int& dau2_flav,
                                              const bool& svfit_conv,
                                              const bool& hh_kinfit_conv,
                                              const float& b_1_hhbtag,
                                              const float& b_2_hhbtag,
                                              const bool& cut_pass,
                                              const float& met_et,
                                              const float& met_phi,
                                              const float& dau1_dxy,
                                              const float& dau2_dxy,
                                              const float& dau1_dz,
                                              const float& dau2_dz,
                                              const float& dau1_iso,
                                              const float& dau2_iso,
                                              const float& DeepMET_ResponseTune_px,
                                              const float& DeepMET_ResponseTune_py,
                                              const float& DeepMET_ResolutionTune_px,
                                              const float& DeepMET_ResolutionTune_py,
                                              const float& met_cov_00,
                                              const float& met_cov_01,
                                              const float& met_cov_11,
                                              const float& bjet1_pnet_bb,
                                              const float& bjet1_pnet_cc,
                                              const float& bjet1_pnet_b,
                                              const float& bjet1_pnet_c,
                                              const float& bjet1_pnet_g,
                                              const float& bjet1_pnet_uds,
                                              const float& bjet1_pnet_pu,
                                              const float& bjet1_pnet_undef,
                                              const float& bjet1_bID,
                                              const float& bjet1_cID,
                                              const float& bjet2_pnet_bb,
                                              const float& bjet2_pnet_cc,
                                              const float& bjet2_pnet_b,
                                              const float& bjet2_pnet_c,
                                              const float& bjet2_pnet_g,
                                              const float& bjet2_pnet_uds,
                                              const float& bjet2_pnet_pu,
                                              const float& bjet2_pnet_undef,
                                              const float& bjet2_bID,
                                              const float& bjet2_cID) {
    /* Processes (requested) features for an event and returns a map of features->values */
    
    std::map<std::string, float> feats = _feat_comp->process(b_1,
                                                             b_2, 
                                                             l_1,
                                                             l_2,
                                                             met,
                                                             svfit,
                                                             vbf_1,
                                                             vbf_2,
                                                             Nu_1,
                                                             Nu_2,
                                                             HHKin_mass_rawass,
                                                             is_boosted,
                                                             pairType,
                                                             year,
                                                             n_vbf,
                                                             dau1_flav,
                                                             dau2_flav,
                                                             svfit_conv,
                                                             hh_kinfit_conv,
                                                             met_et,
                                                             met_phi);
    bool use_vbf = n_vbf >= 2;
    
    // Non-comp extra HL
    if (EvtProc::_feat_check("HHKin_mass_raw_chi2")) feats["HHKin_mass_raw_chi2"] = hh_kinfit_conv ? hh_kinfit_chi2 : -1;
    if (EvtProc::_feat_check("mt2"))            feats["mt2"]            = mt2;
    if (EvtProc::_feat_check("res_mass"))       feats["res_mass"]       = res_mass;
    if (EvtProc::_feat_check("spin"))           feats["spin"]           = spin;
    if (EvtProc::_feat_check("cut_pass"))       feats["cut_pass"]       = cut_pass;

    // Non-comp extra LL
    if (EvtProc::_feat_check("dau1_E"))    feats["dau1_E"]    = l_1.E();
    if (EvtProc::_feat_check("dau2_E"))    feats["dau2_E"]    = l_2.E();
    if (EvtProc::_feat_check("dau1_pt"))   feats["dau1_pt"]   = l_1.Pt();
    if (EvtProc::_feat_check("dau2_pt"))   feats["dau2_pt"]   = l_2.Pt();
    if (EvtProc::_feat_check("bjet1_E"))    feats["bjet1_E"]    = b_1.E();
    if (EvtProc::_feat_check("bjet2_E"))    feats["bjet2_E"]    = b_2.E();
    if (EvtProc::_feat_check("bjet1_pt"))   feats["bjet1_pt"]   = b_1.Pt();
    if (EvtProc::_feat_check("bjet2_pt"))   feats["bjet2_pt"]   = b_2.Pt();
    if (EvtProc::_feat_check("met_pt"))   feats["met_pt"]   = met.Pt();
    if (EvtProc::_feat_check("vbf_1_pt")) feats["vbf_1_pt"] = use_vbf ? vbf_1.Pt() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_pt")) feats["vbf_2_pt"] = use_vbf ? vbf_2.Pt() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_1_E"))  feats["vbf_1_E"]  = use_vbf ? vbf_1.E()  : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_E"))  feats["vbf_2_E"]  = use_vbf ? vbf_2.E()  : std::nanf("1");

    //3-vectors
    if (EvtProc::_feat_check("dau1_px"))   feats["dau1_px"]   = l_1.Px();
    if (EvtProc::_feat_check("dau1_py"))   feats["dau1_py"]   = l_1.Py();
    if (EvtProc::_feat_check("dau1_pz"))   feats["dau1_pz"]   = l_1.Pz();
    if (EvtProc::_feat_check("dau2_px"))   feats["dau2_px"]   = l_2.Px();
    if (EvtProc::_feat_check("dau2_py"))   feats["dau2_py"]   = l_2.Py();
    if (EvtProc::_feat_check("dau2_pz"))   feats["dau2_pz"]   = l_2.Pz();
    if (EvtProc::_feat_check("bjet1_px"))   feats["bjet1_px"]   = b_1.Px();
    if (EvtProc::_feat_check("bjet1_py"))   feats["bjet1_py"]   = b_1.Py();
    if (EvtProc::_feat_check("bjet1_pz"))   feats["bjet1_pz"]   = b_1.Pz();
    if (EvtProc::_feat_check("bjet2_px"))   feats["bjet2_px"]   = b_2.Px();
    if (EvtProc::_feat_check("bjet2_py"))   feats["bjet2_py"]   = b_2.Py();
    if (EvtProc::_feat_check("bjet2_pz"))   feats["bjet2_pz"]   = b_2.Pz();
    if (EvtProc::_feat_check("vbf_1_px")) feats["vbf_1_px"] = use_vbf ? vbf_1.Px() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_1_py")) feats["vbf_1_py"] = use_vbf ? vbf_1.Py() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_1_pz")) feats["vbf_1_pz"] = use_vbf ? vbf_1.Pz() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_px")) feats["vbf_2_px"] = use_vbf ? vbf_2.Px() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_py")) feats["vbf_2_py"] = use_vbf ? vbf_2.Py() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_pz")) feats["vbf_2_pz"] = use_vbf ? vbf_2.Pz() : std::nanf("1");
    if (EvtProc::_feat_check("met_px"))   feats["met_px"]   = met.Px();
    if (EvtProc::_feat_check("met_py"))   feats["met_py"]   = met.Py();
    if (EvtProc::_feat_check("met_E"))    feats["met_E"]    = met.E();

    if (EvtProc::_feat_check("dau1_dxy")) feats["dau1_dxy"] =  dau1_dxy;
    if (EvtProc::_feat_check("dau2_dxy")) feats["dau2_dxy"] =  dau2_dxy;
    if (EvtProc::_feat_check("dau1_dz")) feats["dau1_dz"] =  dau1_dz;
    if (EvtProc::_feat_check("dau2_dz")) feats["dau2_dz"] =  dau2_dz;
    if (EvtProc::_feat_check("dau1_iso")) feats["dau1_iso"] =  dau1_iso;
    if (EvtProc::_feat_check("dau2_iso")) feats["dau2_iso"] =  dau2_iso;

    if (EvtProc::_feat_check("bjet1_pnet_bb")) feats["bjet1_pnet_bb"] =  bjet1_pnet_bb;
    if (EvtProc::_feat_check("bjet1_pnet_cc")) feats["bjet1_pnet_cc"] =  bjet1_pnet_cc;
    if (EvtProc::_feat_check("bjet1_pnet_b")) feats["bjet1_pnet_b"] =  bjet1_pnet_b;
    if (EvtProc::_feat_check("bjet1_pnet_c")) feats["bjet1_pnet_c"] =  bjet1_pnet_c;
    if (EvtProc::_feat_check("bjet1_pnet_g")) feats["bjet1_pnet_g"] =  bjet1_pnet_g;
    if (EvtProc::_feat_check("bjet1_pnet_uds")) feats["bjet1_pnet_uds"] =  bjet1_pnet_uds;
    if (EvtProc::_feat_check("bjet1_pnet_pu")) feats["bjet1_pnet_pu"] =  bjet1_pnet_pu;
    if (EvtProc::_feat_check("bjet1_pnet_undef")) feats["bjet1_pnet_undef"] =  bjet1_pnet_undef;
    if (EvtProc::_feat_check("bjet2_pnet_bb")) feats["bjet2_pnet_bb"] =  bjet2_pnet_bb;
    if (EvtProc::_feat_check("bjet2_pnet_cc")) feats["bjet2_pnet_cc"] =  bjet2_pnet_cc;
    if (EvtProc::_feat_check("bjet2_pnet_b")) feats["bjet2_pnet_b"] =  bjet2_pnet_b;
    if (EvtProc::_feat_check("bjet2_pnet_c")) feats["bjet2_pnet_c"] =  bjet2_pnet_c;
    if (EvtProc::_feat_check("bjet2_pnet_g")) feats["bjet2_pnet_g"] =  bjet2_pnet_g;
    if (EvtProc::_feat_check("bjet2_pnet_uds")) feats["bjet2_pnet_uds"] =  bjet2_pnet_uds;
    if (EvtProc::_feat_check("bjet2_pnet_pu")) feats["bjet2_pnet_pu"] =  bjet2_pnet_pu;
    if (EvtProc::_feat_check("bjet2_pnet_undef")) feats["bjet2_pnet_undef"] =  bjet2_pnet_undef;
    if (EvtProc::_feat_check("bjet1_HHbtag")) feats["bjet1_HHbtag"] =  b_1_hhbtag;
    if (EvtProc::_feat_check("bjet1_bID_deepFlavor")) feats["bjet1_bID_deepFlavor"] = bjet1_bID;
    if (EvtProc::_feat_check("bjet1_cID_deepFlavor")) feats["bjet1_cID_deepFlavor"] = bjet1_cID;
    if (EvtProc::_feat_check("bjet2_HHbtag")) feats["bjet2_HHbtag"] =  b_2_hhbtag;
    if (EvtProc::_feat_check("bjet2_bID_deepFlavor")) feats["bjet2_bID_deepFlavor"] = bjet2_bID;
    if (EvtProc::_feat_check("bjet2_cID_deepFlavor")) feats["bjet2_cID_deepFlavor"] = bjet2_cID;
    if (EvtProc::_feat_check("dmet_resp_px")) feats["dmet_resp_px"] = DeepMET_ResponseTune_px;
    if (EvtProc::_feat_check("dmet_resp_py")) feats["dmet_resp_py"] = DeepMET_ResponseTune_py;
    if (EvtProc::_feat_check("dmet_reso_px")) feats["dmet_reso_px"] = DeepMET_ResolutionTune_px;
    if (EvtProc::_feat_check("dmet_reso_py")) feats["dmet_reso_py"] = DeepMET_ResolutionTune_py;
    if (EvtProc::_feat_check("met_cov00")) feats["met_cov00"] = met_cov_00;
    if (EvtProc::_feat_check("met_cov01")) feats["met_cov01"] = met_cov_01;
    if (EvtProc::_feat_check("met_cov11")) feats["met_cov11"] = met_cov_11;

    // Htautau inputs
    if (EvtProc::_feat_check("pairType"))          feats["pairType"]          = pairType;
    if (EvtProc::_feat_check("dau1_decayMode"))    feats["dau1_decayMode"]    = dau1_decayMode;
    if (EvtProc::_feat_check("dau2_decayMode"))    feats["dau2_decayMode"]    = dau2_decayMode;

    return feats;
}

inline bool EvtProc::_feat_check(std::string feat) {return (_all ? true : std::find(_requested.begin(), _requested.end(), feat) != _requested.end());}

std::vector<float> EvtProc::process_as_vec(const LorentzVector& b_1,
                                              const LorentzVector& b_2,
                                              const LorentzVector& l_1,
                                              const LorentzVector& l_2,
                                              const LorentzVector& met,
                                              const LorentzVector& svfit,
                                              const LorentzVector& vbf_1,
                                              const LorentzVector& vbf_2,
                                              const LorentzVector& Nu_1,
                                              const LorentzVector& Nu_2,
                                              const float& HHKin_mass_rawass,
                                              const float& hh_kinfit_chi2,
                                              const float& mt2,
                                              const bool&  is_boosted,
                                              Year year,
                                              const float& res_mass,
                                              Spin spin,
                                              const int& n_vbf,
                                              const int& pairType,
                                              const int& dau1_decayMode,
                                              const int& dau2_decayMode,
                                              const int& dau1_flav,
                                              const int& dau2_flav,
                                              const bool& svfit_conv,
                                              const bool& hh_kinfit_conv,
                                              const float& b_1_hhbtag,
                                              const float& b_2_hhbtag,
                                              const bool& cut_pass,
                                              const float& met_et,
                                              const float& met_phi,
                                              const float& dau1_dxy,
                                              const float& dau2_dxy,
                                              const float& dau1_dz,
                                              const float& dau2_dz,
                                              const float& dau1_iso,
                                              const float& dau2_iso,
                                              const float& DeepMET_ResponseTune_px,
                                              const float& DeepMET_ResponseTune_py,
                                              const float& DeepMET_ResolutionTune_px,
                                              const float& DeepMET_ResolutionTune_py,
                                              const float& met_cov_00,
                                              const float& met_cov_01,
                                              const float& met_cov_11,
                                              const float& bjet1_pnet_bb,
                                              const float& bjet1_pnet_cc,
                                              const float& bjet1_pnet_b,
                                              const float& bjet1_pnet_c,
                                              const float& bjet1_pnet_g,
                                              const float& bjet1_pnet_uds,
                                              const float& bjet1_pnet_pu,
                                              const float& bjet1_pnet_undef,
                                              const float& bjet1_bID,
                                              const float& bjet1_cID,
                                              const float& bjet2_pnet_bb,
                                              const float& bjet2_pnet_cc,
                                              const float& bjet2_pnet_b,
                                              const float& bjet2_pnet_c,
                                              const float& bjet2_pnet_g,
                                              const float& bjet2_pnet_uds,
                                              const float& bjet2_pnet_pu,
                                              const float& bjet2_pnet_undef,
                                              const float& bjet2_bID,
                                              const float& bjet2_cID){
    /* Calls  EvtProc::process and processes result into a vector */

    std::map<std::string, float> feats = EvtProc::process(b_1,
                                              b_2,
                                              l_1,
                                              l_2,
                                              met,
                                              svfit,
                                              vbf_1,
                                              vbf_2,
                                              Nu_1,
                                              Nu_2,
                                              HHKin_mass_rawass,
                                              hh_kinfit_chi2,
                                              mt2,
                                              is_boosted,
                                              year,
                                              res_mass,
                                              spin,
                                              n_vbf,
                                              pairType,
                                              dau1_decayMode,
                                              dau2_decayMode,
                                              dau1_flav,
                                              dau2_flav,
                                              svfit_conv,
                                              hh_kinfit_conv,
                                              b_1_hhbtag,
                                              b_2_hhbtag,
                                              cut_pass,
                                              met_et,
                                              met_phi,
                                              dau1_dxy,
                                              dau2_dxy,
                                              dau1_dz,
                                              dau2_dz,
                                              dau1_iso,
                                              dau2_iso,
                                              DeepMET_ResponseTune_px,
                                              DeepMET_ResponseTune_py,
                                              DeepMET_ResolutionTune_px,
                                              DeepMET_ResolutionTune_py,
                                              met_cov_00,
                                              met_cov_01,
                                              met_cov_11,
                                              bjet1_pnet_bb,
                                              bjet1_pnet_cc,
                                              bjet1_pnet_b,
                                              bjet1_pnet_c,
                                              bjet1_pnet_g,
                                              bjet1_pnet_uds,
                                              bjet1_pnet_pu,
                                              bjet1_pnet_undef,
                                              bjet1_bID,
                                              bjet1_cID,
                                              bjet2_pnet_bb,
                                              bjet2_pnet_cc,
                                              bjet2_pnet_b,
                                              bjet2_pnet_c,
                                              bjet2_pnet_g,
                                              bjet2_pnet_uds,
                                              bjet2_pnet_pu,
                                              bjet2_pnet_undef,
                                              bjet2_bID,
                                              bjet2_cID);
    std::vector<float> vec(feats.size());
    int i = 0;
    if (_all) {
        for (auto const& f : feats) {
            vec[i] = f.second;
            i++;
        }
    } else {
        for (auto const& f : _requested) {
            vec[i] = feats[f];
            i++;
        }
    }
    return vec;
}

void EvtProc::process_to_vec(std::vector<std::unique_ptr<float>>& feats,
                                              const LorentzVector& b_1,
                                              const LorentzVector& b_2,
                                              const LorentzVector& l_1,
                                              const LorentzVector& l_2,
                                              const LorentzVector& met,
                                              const LorentzVector& svfit,
                                              const LorentzVector& vbf_1,
                                              const LorentzVector& vbf_2,
                                              const LorentzVector& Nu_1,
                                              const LorentzVector& Nu_2,
                                              const float& HHKin_mass_rawass,
                                              const float& hh_kinfit_chi2,
                                              const float& mt2,
                                              const bool&  is_boosted,
                                              Year year,
                                              const float& res_mass,
                                              Spin spin,
                                              const int& n_vbf,
                                              const int& pairType,
                                              const int& dau1_decayMode,
                                              const int& dau2_decayMode,
                                              const int& dau1_flav,
                                              const int& dau2_flav,
                                              const bool& svfit_conv,
                                              const bool& hh_kinfit_conv,
                                              const float& b_1_hhbtag,
                                              const float& b_2_hhbtag,
                                              const bool& cut_pass,
                                              const float& met_et,
                                              const float& met_phi,
                                              const float& dau1_dxy,
                                              const float& dau2_dxy,
                                              const float& dau1_dz,
                                              const float& dau2_dz,
                                              const float& dau1_iso,
                                              const float& dau2_iso,
                                              const float& DeepMET_ResponseTune_px,
                                              const float& DeepMET_ResponseTune_py,
                                              const float& DeepMET_ResolutionTune_px,
                                              const float& DeepMET_ResolutionTune_py,
                                              const float& met_cov_00,
                                              const float& met_cov_01,
                                              const float& met_cov_11,
                                              const float& bjet1_pnet_bb,
                                              const float& bjet1_pnet_cc,
                                              const float& bjet1_pnet_b,
                                              const float& bjet1_pnet_c,
                                              const float& bjet1_pnet_g,
                                              const float& bjet1_pnet_uds,
                                              const float& bjet1_pnet_pu,
                                              const float& bjet1_pnet_undef,
                                              const float& bjet1_bID,
                                              const float& bjet1_cID,
                                              const float& bjet2_pnet_bb,
                                              const float& bjet2_pnet_cc,
                                              const float& bjet2_pnet_b,
                                              const float& bjet2_pnet_c,
                                              const float& bjet2_pnet_g,
                                              const float& bjet2_pnet_uds,
                                              const float& bjet2_pnet_pu,
                                              const float& bjet2_pnet_undef,
                                              const float& bjet2_bID,
                                              const float& bjet2_cID) {
    /* Calls  EvtProc::process and processes result into a supplied vector */

    std::map<std::string, float> feat_vals = EvtProc::process(b_1,
                                              b_2,
                                              l_1,
                                              l_2,
                                              met,
                                              svfit,
                                              vbf_1,
                                              vbf_2,
                                              Nu_1,
                                              Nu_2,
                                              HHKin_mass_rawass,
                                              hh_kinfit_chi2,
                                              mt2,
                                              is_boosted,
                                              year,
                                              res_mass,
                                              spin,
                                              n_vbf,
                                              pairType,
                                              dau1_decayMode,
                                              dau2_decayMode,
                                              dau1_flav,
                                              dau2_flav,
                                              svfit_conv,
                                              hh_kinfit_conv,
                                              b_1_hhbtag,
                                              b_2_hhbtag,
                                              cut_pass,
                                              met_et,
                                              met_phi,
                                              dau1_dxy,
                                              dau2_dxy,
                                              dau1_dz,
                                              dau2_dz,
                                              dau1_iso,
                                              dau2_iso,
                                              DeepMET_ResponseTune_px,
                                              DeepMET_ResponseTune_py,
                                              DeepMET_ResolutionTune_px,
                                              DeepMET_ResolutionTune_py,
                                              met_cov_00,
                                              met_cov_01,
                                              met_cov_11,
                                              bjet1_pnet_bb,
                                              bjet1_pnet_cc,
                                              bjet1_pnet_b,
                                              bjet1_pnet_c,
                                              bjet1_pnet_g,
                                              bjet1_pnet_uds,
                                              bjet1_pnet_pu,
                                              bjet1_pnet_undef,
                                              bjet1_bID,
                                              bjet1_cID,
                                              bjet2_pnet_bb,
                                              bjet2_pnet_cc,
                                              bjet2_pnet_b,
                                              bjet2_pnet_c,
                                              bjet2_pnet_g,
                                              bjet2_pnet_uds,
                                              bjet2_pnet_pu,
                                              bjet2_pnet_undef,
                                              bjet2_bID,
                                              bjet2_cID);
    if (feat_vals.size() != feats.size()) throw std::length_error("Length of computed map (" + std::to_string(feat_vals.size()) + ") does not match length of \
                                                                   vector to fill (" + std::to_string(feats.size()) + ")\n");
    int i = 0;
    if (_all) {
        for (auto const& f : feat_vals) {
            *(feats[i]) = f.second;
            i++;
        }
    } else {
        for (auto const& f : _requested) {
            *(feats[i]) = feat_vals[f];
            i++;
        }
    }
    
}

std::vector<std::string> EvtProc::get_feats() {
    /* Returns list of features that will be computed in general operation */

    std::map<std::string, float> feats = EvtProc::process(LorentzVector(),
                                              LorentzVector(),
                                              LorentzVector(),
                                              LorentzVector(),
                                              LorentzVector(),
                                              LorentzVector(),
                                              LorentzVector(),
                                              LorentzVector(),
                                              LorentzVector(),
                                              LorentzVector(),
                                              0,
                                              0,
                                              0,
                                              true,
                                              Year(y18),
                                              0,
                                              Spin(nonres),
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              true,
                                              true,
                                              0,
                                              0,
                                              true,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0,
                                              0);
    std::vector<std::string> names;
    if (_all) {
        for (auto const& f : feats) names.push_back(f.first);
    } else {
        for (auto const& f : _requested) {
            if (feats.find(f) == feats.end()) throw std::invalid_argument("Requested feature " + f + " not compued\n");
            names.push_back(f);
        }
    }
    return names;
}