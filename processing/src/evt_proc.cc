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
                                              const float& hh_kinfit_mass,
                                              const float& hh_kinfit_chi2,
                                              const float& mt2,
                                              const bool&  is_boosted,
                                              const float& b_1_csv,
                                              const float& b_2_csv,
                                              Channel channel,
                                              Year year,
                                              const float& res_mass,
                                              Spin spin,
                                              const float& klambda,
                                              const int& n_vbf,
                                              const bool& svfit_conv,
                                              const bool& hh_kinfit_conv,
                                              const float& b_1_hhbtag,
										      const float& b_2_hhbtag,
										      const float& vbf_1_hhbtag,
										      const float& vbf_2_hhbtag,
                                              const float& b_1_cvsl,
										      const float& b_2_cvsl,
										      const float& vbf_1_cvsl,
										      const float& vbf_2_cvsl,
										      const float& b_1_cvsb,
										      const float& b_2_cvsb,
										      const float& vbf_1_cvsb,
										      const float& vbf_2_cvsb,
                                              const float& cv,
										      const float& c2v,
										      const float& c3) {
    /* Processes (requested) features for an event and returns a map of features->values */
    
    std::map<std::string, float> feats = _feat_comp->process(b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, hh_kinfit_mass, is_boosted,
                                                             b_1_csv, b_2_csv, channel, year, n_vbf, svfit_conv, hh_kinfit_conv);
    bool use_vbf = n_vbf >= 2;
    
    // Non-comp extra HL
    if (EvtProc::_feat_check("hh_kinfit_chi2")) feats["hh_kinfit_chi2"] = hh_kinfit_conv ? hh_kinfit_chi2 : std::nanf("1");
    if (EvtProc::_feat_check("mt2"))            feats["mt2"]            = mt2;
    if (EvtProc::_feat_check("mt_tot"))         feats["mt_tot"]         = mt_tot; 
    if (EvtProc::_feat_check("res_mass"))       feats["res_mass"]       = res_mass;
    if (EvtProc::_feat_check("spin"))           feats["spin"]           = spin;
    if (EvtProc::_feat_check("klambda"))        feats["klambda"]        = klambda;
    if (EvtProc::_feat_check("cv"))             feats["cv"]             = use_vbf ? cv  : std::nanf("1");
    if (EvtProc::_feat_check("c2v"))            feats["c2v"]            = use_vbf ? c2v : std::nanf("1");
    if (EvtProc::_feat_check("c3"))             feats["c3"]             = use_vbf ? c3  : std::nanf("1");
    if (EvtProc::_feat_check("b_1_hhbtag"))     feats["b_1_hhbtag"]     = b_1_hhbtag;
    if (EvtProc::_feat_check("b_2_hhbtag"))     feats["b_2_hhbtag"]     = b_2_hhbtag;
    if (EvtProc::_feat_check("vbf_1_hhbtag"))   feats["vbf_1_hhbtag"]   = use_vbf ? vbf_1_hhbtag : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_hhbtag"))   feats["vbf_2_hhbtag"]   = use_vbf ? vbf_2_hhbtag : std::nanf("1");
    if (EvtProc::_feat_check("b_1_cvsl"))       feats["b_1_cvsl"]       = b_1_cvsl;
    if (EvtProc::_feat_check("b_2_cvsl"))       feats["b_2_cvsl"]       = b_2_cvsl;
    if (EvtProc::_feat_check("vbf_1_cvsl"))     feats["vbf_1_cvsl"]     = use_vbf ? vbf_1_cvsl : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_cvsl"))     feats["vbf_2_cvsl"]     = use_vbf ? vbf_2_cvsl : std::nanf("1");
    if (EvtProc::_feat_check("b_1_cvsb"))       feats["b_1_cvsb"]       = b_1_cvsb;
    if (EvtProc::_feat_check("b_2_cvsb"))       feats["b_2_cvsb"]       = b_2_cvsb;
    if (EvtProc::_feat_check("vbf_1_cvsb"))     feats["vbf_1_cvsb"]     = use_vbf ? vbf_1_cvsb : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_cvsb"))     feats["vbf_2_cvsb"]     = use_vbf ? vbf_2_cvsb : std::nanf("1");

    // Non-comp extra LL
    if (EvtProc::_feat_check("l_1_E"))    feats["l_1_E"]    = l_1.E();
    if (EvtProc::_feat_check("l_2_E"))    feats["l_2_E"]    = l_2.E();
    if (EvtProc::_feat_check("l_1_pT"))   feats["l_1_pT"]   = l_1.Pt();
    if (EvtProc::_feat_check("l_2_pT"))   feats["l_2_pT"]   = l_2.Pt();
    if (EvtProc::_feat_check("b_1_E"))    feats["b_1_E"]    = b_1.E();
    if (EvtProc::_feat_check("b_2_E"))    feats["b_2_E"]    = b_2.E();
    if (EvtProc::_feat_check("b_1_pT"))   feats["b_1_pT"]   = b_1.Pt();
    if (EvtProc::_feat_check("b_2_pT"))   feats["b_2_pT"]   = b_2.Pt();
    if (EvtProc::_feat_check("met_pT"))   feats["met_pT"]   = met.Pt();
    if (EvtProc::_feat_check("vbf_1_pT")) feats["vbf_1_pT"] = use_vbf ? vbf_1.Pt() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_pT")) feats["vbf_2_pT"] = use_vbf ? vbf_2.Pt() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_1_E"))  feats["vbf_1_E"]  = use_vbf ? vbf_1.E()  : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_E"))  feats["vbf_2_E"]  = use_vbf ? vbf_2.E()  : std::nanf("1");

    //3-vectors
    if (EvtProc::_feat_check("l_1_px"))   feats["l_1_px"]   = l_1.Px();
    if (EvtProc::_feat_check("l_1_py"))   feats["l_1_py"]   = l_1.Py();
    if (EvtProc::_feat_check("l_1_pz"))   feats["l_1_pz"]   = l_1.Pz();
    if (EvtProc::_feat_check("l_2_px"))   feats["l_2_px"]   = l_2.Px();
    if (EvtProc::_feat_check("l_2_py"))   feats["l_2_py"]   = l_2.Py();
    if (EvtProc::_feat_check("l_2_pz"))   feats["l_2_pz"]   = l_2.Pz();
    if (EvtProc::_feat_check("b_1_px"))   feats["b_1_px"]   = b_1.Px();
    if (EvtProc::_feat_check("b_1_py"))   feats["b_1_py"]   = b_1.Py();
    if (EvtProc::_feat_check("b_1_pz"))   feats["b_1_pz"]   = b_1.Pz();
    if (EvtProc::_feat_check("b_2_px"))   feats["b_2_px"]   = b_2.Px();
    if (EvtProc::_feat_check("b_2_py"))   feats["b_2_py"]   = b_2.Py();
    if (EvtProc::_feat_check("b_2_pz"))   feats["b_2_pz"]   = b_2.Pz();
    if (EvtProc::_feat_check("vbf_1_px")) feats["vbf_1_px"] = use_vbf ? vbf_1.Px() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_1_py")) feats["vbf_1_py"] = use_vbf ? vbf_1.Py() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_1_pz")) feats["vbf_1_pz"] = use_vbf ? vbf_1.Pz() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_px")) feats["vbf_2_px"] = use_vbf ? vbf_2.Px() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_py")) feats["vbf_2_py"] = use_vbf ? vbf_2.Py() : std::nanf("1");
    if (EvtProc::_feat_check("vbf_2_pz")) feats["vbf_2_pz"] = use_vbf ? vbf_2.Pz() : std::nanf("1");
    if (EvtProc::_feat_check("met_px"))   feats["met_px"]   = met.Px();
    if (EvtProc::_feat_check("met_py"))   feats["met_py"]   = met.Py();
    if (EvtProc::_feat_check("met_E"))    feats["met_E"]    = met.E();
    
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
                                           const float& hh_kinfit_mass,
                                           const float& hh_kinfit_chi2,
                                           const float& mt2,
                                           const bool& is_boosted,
                                           const float& b_1_csv,
                                           const float& b_2_csv,
                                           Channel channel,
                                           Year year,
                                           const float& res_mass,
                                           Spin spin,
                                           const float& klambda,
                                           const int& n_vbf,
                                           const bool& svfit_conv,
                                           const bool& hh_kinfit_conv,
                                           const float& b_1_hhbtag,
										   const float& b_2_hhbtag,
										   const float& vbf_1_hhbtag,
										   const float& vbf_2_hhbtag,
                                           const float& b_1_cvsl,
										   const float& b_2_cvsl,
										   const float& vbf_1_cvsl,
										   const float& vbf_2_cvsl,
										   const float& b_1_cvsb,
										   const float& b_2_cvsb,
										   const float& vbf_1_cvsb,
										   const float& vbf_2_cvsb,
                                           const float& cv,
										   const float& c2v,
										   const float& c3) {
    /* Calls  EvtProc::process and processes result into a vector */

    std::map<std::string, float> feats = EvtProc::process(b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, hh_kinfit_mass, hh_kinfit_chi2, mt2, is_boosted,
                                                          b_1_csv, b_2_csv, channel, year, res_mass, spin, klambda, n_vbf, svfit_conv, hh_kinfit_conv,
                                                          b_1_hhbtag, b_2_hhbtag, vbf_1_hhbtag, vbf_2_hhbtag,
                                                          b_1_cvsl, b_2_cvsl, vbf_1_cvsl, vbf_2_cvsl, b_1_cvsb, b_2_cvsb, vbf_1_cvsb, vbf_2_cvsb,
                                                          cv, c2v, c3);
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
                             const float& hh_kinfit_mass,
                             const float& hh_kinfit_chi2,
                             const float& mt2,
                             const bool& is_boosted,
                             const float& b_1_csv,
                             const float& b_2_csv,
                             Channel channel,
                             Year year,
                             const float& res_mass,
                             Spin spin,
                             const float& klambda,
                             const int& n_vbf,
                             const bool& svfit_conv,
                             const bool& hh_kinfit_conv,
                             const float& b_1_hhbtag,
							 const float& b_2_hhbtag,
							 const float& vbf_1_hhbtag,
							 const float& vbf_2_hhbtag,
                             const float& b_1_cvsl,
							 const float& b_2_cvsl,
							 const float& vbf_1_cvsl,
							 const float& vbf_2_cvsl,
							 const float& b_1_cvsb,
							 const float& b_2_cvsb,
							 const float& vbf_1_cvsb,
							 const float& vbf_2_cvsb,
                             const float& cv,
							 const float& c2v,
						     const float& c3) {
    /* Calls  EvtProc::process and processes result into a supplied vector */

    std::map<std::string, float> feat_vals = EvtProc::process(b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, hh_kinfit_mass, hh_kinfit_chi2, mt2,
                                                              is_boosted, b_1_csv, b_2_csv, channel, year, res_mass, spin, klambda, n_vbf, svfit_conv,
                                                              hh_kinfit_conv, b_1_hhbtag, b_2_hhbtag, vbf_1_hhbtag, vbf_2_hhbtag,
                                                              b_1_cvsl, b_2_cvsl, vbf_1_cvsl, vbf_2_cvsl, b_1_cvsb, b_2_cvsb, vbf_1_cvsb, vbf_2_cvsb,
                                                              cv, c2v, c3);
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

    std::map<std::string, float> feats = EvtProc::process(LorentzVector(), LorentzVector(), LorentzVector(), LorentzVector(), LorentzVector(), LorentzVector(),
                                                          LorentzVector(), LorentzVector(), 0, 0, 0, 0, 0, 0, 0, 0, false, 0, 0, 0, 0, Channel(tauTau),
                                                          Year(y16), 0, Spin(nonres), 0, 2, true, true, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
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