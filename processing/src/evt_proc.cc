#include "cms_hh_proc_interface/processing/interface/evt_proc.hh"

EvtProc::EvtProc(bool return_all, std::set<std::string> requested, bool use_deep_csv) {
    _all = return_all;
    _requested = requested;
    _feat_comp = new FeatComp(return_all, requested, use_deep_csv);
}

EvtProc::~EvtProc() {}

std::map<std::string, float> EvtProc::process(const LorentzVector& b_1,
                                              const LorentzVector& b_2,
                                              const LorentzVector& l_1,
                                              const LorentzVector& l_2,
                                              const LorentzVector& met,
                                              const LorentzVector& svfit,
                                              const float& hh_kinfit_mass,
                                              const float& hh_kinfit_chi2,
                                              const float& mt2,
                                              const float& mt_tot,
                                              const bool& is_boosted,
                                              const float& b_1_csv,
                                              const float& b_2_csv,
                                              const float& b_1_deepcsv,
                                              const float& b_2_deepcsv,
                                              Channel channel,
                                              const float& res_mass) {
    /* Processes (requested) features for an event and returns a map of features->values (in order of requested features) */

    std::map<std::string, float> feats = _feat_comp->process(b_1, b_2, l_1, l_2, met, svfit, hh_kinfit_mass, is_boosted,
                                                             b_1_csv, b_2_csv, b_1_deepcsv, b_2_deepcsv, channel);
    if (EvtProc::_feat_check("hh_kinfit_chi2")) feats["hh_kinfit_chi2"] = hh_kinfit_chi2;
    if (EvtProc::_feat_check("mt2"))            feats["mt2"]            = mt2;
    if (EvtProc::_feat_check("mt_tot"))         feats["mt_tot"]         = mt_tot; 
    if (EvtProc::_feat_check("res_mass"))       feats["res_mass"]       = res_mass;

    return _all ? feats : EvtProc::_sort_feats(feats);
}

inline bool EvtProc::_feat_check(std::string feat) {return (_all ? true : _requested.find(feat) != _requested.end());}

std::map<std::string, float> EvtProc::_sort_feats(std::map<std::string, float> feats) {
    std::map<std::string, float> sf;
    for (auto const& f : _requested) sf[f] = feats[f];
    return sf;
}
