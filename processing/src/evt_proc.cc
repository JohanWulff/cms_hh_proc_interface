#include "cms_hh_proc_interface/processing/interface/evt_proc.hh"

EvtProc::EvtProc(bool return_all, std::vector<std::string> requested, bool use_deep_csv) {
    _all = return_all;
    _requested = requested;
    _feat_comp = new FeatComp(return_all, _requested, use_deep_csv);
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
                                              const float& mt_tot,
                                              const float& p_zetavisible,
                                              const float& p_zeta,
                                              const float& top_1_mass,
                                              const float& top_2_mass,
                                              const float& l_1_mt,
                                              const float& l_2_mt,
                                              const bool&  is_boosted,
                                              const float& b_1_csv,
                                              const float& b_2_csv,
                                              const float& b_1_deepcsv,
                                              const float& b_2_deepcsv,
                                              Channel channel,
                                              Year year,
                                              const float& res_mass,
                                              Spin spin,
                                              const float& klambda) {
    /* Processes (requested) features for an event and returns a map of features->values */
    
    std::map<std::string, float> feats = _feat_comp->process(b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, hh_kinfit_mass, is_boosted,
                                                             b_1_csv, b_2_csv, b_1_deepcsv, b_2_deepcsv, channel, year);
    // Non-comp extra HL
    if (EvtProc::_feat_check("hh_kinfit_chi2")) feats["hh_kinfit_chi2"] = hh_kinfit_chi2 >= 0 ? hh_kinfit_chi2 : std::nanf("1");
    if (EvtProc::_feat_check("mt2"))            feats["mt2"]            = mt2;
    if (EvtProc::_feat_check("mt_tot"))         feats["mt_tot"]         = mt_tot; 
    if (EvtProc::_feat_check("res_mass"))       feats["res_mass"]       = res_mass;
    if (EvtProc::_feat_check("spin"))           feats["spin"]           = spin;
    if (EvtProc::_feat_check("klambda"))        feats["klambda"]        = klambda;
    if (EvtProc::_feat_check("p_zetavisible"))  feats["p_zetavisible"]  = p_zetavisible;
    if (EvtProc::_feat_check("p_zeta"))         feats["p_zeta"]         = p_zeta;
    if (EvtProc::_feat_check("top_1_mass"))     feats["top_1_mass"]     = top_1_mass;
    if (EvtProc::_feat_check("top_2_mass"))     feats["top_2_mass"]     = top_2_mass;
    
    // Non-comp extra LL
    if (EvtProc::_feat_check("l_1_mt"))   feats["l_1_mt"]   = l_1_mt;
    if (EvtProc::_feat_check("l_2_mt"))   feats["l_2_mt"]   = l_2_mt;
    if (EvtProc::_feat_check("l_1_E"))    feats["l_1_E"]    = l_1.E();
    if (EvtProc::_feat_check("l_2_E"))    feats["l_2_E"]    = l_2.E();
    if (EvtProc::_feat_check("l_1_pT"))   feats["l_1_pT"]   = l_1.Pt();
    if (EvtProc::_feat_check("l_2_pT"))   feats["l_2_pT"]   = l_2.Pt();
    if (EvtProc::_feat_check("b_1_E"))    feats["b_1_E"]    = b_1.E();
    if (EvtProc::_feat_check("b_2_E"))    feats["b_2_E"]    = b_2.E();
    if (EvtProc::_feat_check("b_1_pT"))   feats["b_1_pT"]   = b_1.Pt();
    if (EvtProc::_feat_check("b_2_pT"))   feats["b_2_pT"]   = b_2.Pt();
    if (EvtProc::_feat_check("met_pT"))   feats["met_pT"]   = met.Pt();
    if (EvtProc::_feat_check("vbf_1_pT")) feats["vbf_1_pT"] = vbf_1.Pt();
    if (EvtProc::_feat_check("vbf_2_pT")) feats["vbf_2_pT"] = vbf_2.Pt();
    if (EvtProc::_feat_check("vbf_1_E"))  feats["vbf_1_E"]  = vbf_1.E();
    if (EvtProc::_feat_check("vbf_2_E"))  feats["vbf_2_E"]  = vbf_2.E();

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
                                           const float& mt_tot,
                                           const float& p_zetavisible,
                                           const float& p_zeta,
                                           const float& top_1_mass,
                                           const float& top_2_mass,
                                           const float& l_1_mt,
                                           const float& l_2_mt,
                                           const bool& is_boosted,
                                           const float& b_1_csv,
                                           const float& b_2_csv,
                                           const float& b_1_deepcsv,
                                           const float& b_2_deepcsv,
                                           Channel channel,
                                           Year year,
                                           const float& res_mass,
                                           Spin spin,
                                           const float& klambda) {
    /* Calls  EvtProc::process and processes result into a vector */

    std::map<std::string, float> feats = EvtProc::process(b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, hh_kinfit_mass, hh_kinfit_chi2, mt2, mt_tot, p_zetavisible, p_zeta,
                                                          top_1_mass, top_2_mass, l_1_mt, l_2_mt,
                                                          is_boosted, b_1_csv, b_2_csv, b_1_deepcsv, b_2_deepcsv, channel, year, res_mass, spin, klambda);
    std::cout << "EvtProc::process_as_vec " << feats.size() << _requested.size() << "\n";
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
                             const float& mt_tot,
                             const float& p_zetavisible,
                             const float& p_zeta,
                             const float& top_1_mass,
                             const float& top_2_mass,
                             const float& l_1_mt,
                             const float& l_2_mt,
                             const bool& is_boosted,
                             const float& b_1_csv,
                             const float& b_2_csv,
                             const float& b_1_deepcsv,
                             const float& b_2_deepcsv,
                             Channel channel,
                             Year year,
                             const float& res_mass,
                             Spin spin,
                             const float& klambda) {
    /* Calls  EvtProc::process and processes result into a supplied vector */

    std::map<std::string, float> feat_vals = EvtProc::process(b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, hh_kinfit_mass, hh_kinfit_chi2, mt2, mt_tot, p_zetavisible, p_zeta,
                                                              top_1_mass, top_2_mass, l_1_mt, l_2_mt,
                                                              is_boosted, b_1_csv, b_2_csv, b_1_deepcsv, b_2_deepcsv, channel, year, res_mass, spin, klambda);
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
                                                          LorentzVector(), LorentzVector(), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, false, 0, 0, 0, 0, Channel(tauTau),
                                                          Year(y16), 0, Spin(nonres), 0);
    std::vector<std::string> names;
    for (auto const& f : feats) names.push_back(f.first);
    return names;
}