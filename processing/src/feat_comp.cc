#include "cms_hh_proc_interface/processing/interface/feat_comp.hh"

FeatComp::FeatComp(bool return_all, std::vector<std::string> requested={}, bool use_deep_csv=true, bool verbose=false) {
    _all = return_all;
    _requested = requested;
    _use_deep_csv = use_deep_csv;
    _verbose = verbose;
}

FeatComp::~FeatComp() {}

std::map<std::string, float> FeatComp::process(const TLorentzVector& b_1,
										       const TLorentzVector& b_2, 
										       const TLorentzVector& l_1,
										       const TLorentzVector& l_2,
										       const TLorentzVector& met,
										       const TLorentzVector& svfit,
										       const float& hh_kinfit_m,
                                               const bool& is_boosted,
                                               const float& b_1_csv,
                                               const float& b_2_csv,
                                               const float& b_1_deepcsv,
                                               const float& b_2_deepcsv,
                                               FeatComp::Channel channel) {
    /* Compute HL features from base event*/

    // Extra vectors
    TLorentzVector h_bb(b_1 + b_2);
    TLorentzVector h_tt_vis(l_1 + l_2);
    TLorentzVector h_tt_met(h_tt_vis + met);
    TLorentzVector hh;
    hh.SetXYZM(h_bb.Px()+svfit.Px(), h_bb.Py()+svfit.Py(), h_bb.Pz()+svfit.Pz(), hh_kinfit_m);

    std::map<std::string, float> feats;

    // Categoricals
    feats["boosted"] = is_boosted;
    feats["channel"] = channel;
    FeatComp::_add_jet_flags(b_1_csv, b_2_csv, b_1_deepcsv, b_2_deepcsv, feats);

    // Delta phi
    if (FeatComp::_feat_check("dphi_l1_l2"))      feats["dphi_l1_l2"]      = FeatComp::delta_phi(l_1, l_2);
    if (FeatComp::_feat_check("dphi_b1_b2"))      feats["dphi_b1_b2"]      = FeatComp::delta_phi(b_1, b_2);
    if (FeatComp::_feat_check("dphi_l1_met"))     feats["dphi_l1_met"]     = FeatComp::delta_phi(l_1, met);
    if (FeatComp::_feat_check("dphi_l2_met"))     feats["dphi_l2_met"]     = FeatComp::delta_phi(l_2, met);
    if (FeatComp::_feat_check("dphi_l1l2_met"))   feats["dphi_l1l2_met"]   = FeatComp::delta_phi(h_tt_vis, met);
    if (FeatComp::_feat_check("dphi_sv_met"))     feats["dphi_sv_met"]     = FeatComp::delta_phi(svfit, met);
    if (FeatComp::_feat_check("dphi_httmet_met")) feats["dphi_httmet_met"] = FeatComp::delta_phi(h_tt_met, met);
    if (FeatComp::_feat_check("dphi_httvis_met")) feats["dphi_httvis_met"] = FeatComp::delta_phi(h_tt_vis, met);
    if (FeatComp::_feat_check("dphi_hbb_met"))    feats["dphi_hbb_met"]    = FeatComp::delta_phi(h_bb, met);
    if (FeatComp::_feat_check("dphi_hbb_sv"))     feats["dphi_hbb_sv"]     = FeatComp::delta_phi(h_bb, svfit);
    if (FeatComp::_feat_check("dphi_hbb_httmet")) feats["dphi_hbb_httmet"] = FeatComp::delta_phi(h_bb, h_tt_met);

    // Delta eta
    if (FeatComp::_feat_check("deta_l1_l2"))      feats["deta_l1_l2"]      = FeatComp::delta_eta(l_1, l_2);
    if (FeatComp::_feat_check("deta_b1_b2"))      feats["deta_l1_l2"]      = FeatComp::delta_eta(b_1, b_2);
    if (FeatComp::_feat_check("deta_hbb_sv"))     feats["deta_hbb_sv"]     = FeatComp::delta_eta(h_bb, svfit);
    if (FeatComp::_feat_check("deta_hbb_httmet")) feats["deta_hbb_httmet"] = FeatComp::delta_eta(h_bb, met);

    // Delta R
    if (FeatComp::_feat_check("dR_l1_l2"))      feats["dR_l1_l2"]      = l_1.DeltaR(l_2);
    if (FeatComp::_feat_check("dR_b1_b2"))      feats["dR_b1_b2"]      = b_1.DeltaR(b_2);
    if (FeatComp::_feat_check("dR_hbb_httmet")) feats["dR_hbb_httmet"] = h_bb.DeltaR(h_tt_met);
    if (FeatComp::_feat_check("dR_hbb_sv"))     feats["dR_hbb_sv"]     = h_bb.DeltaR(svfit);
    if (FeatComp::_feat_check("dR_b1_b2_x_h_bb_pT"))     feats["dR_b1_b2_x_h_bb_pT"]     = b_1.DeltaR(b_2)*h_bb.Pt();
    if (FeatComp::_feat_check("dR_l1_l2_x_h_tt_met_pT")) feats["dR_l1_l2_x_h_tt_met_pT"] = l_1.DeltaR(l_2)*h_tt_met.Pt();
    if (FeatComp::_feat_check("dR_l1_l2_x_sv_pT"))       feats["dR_l1_l2_x_sv_pT"]       = l_1.DeltaR(l_2)*svfit.Pt();
    if (FeatComp::_feat_check("dR_b1_b2_boosted_hbb"))     feats["dR_b1_b2_boosted_hbb"]     = FeatComp::delta_r_boosted(b_1, b_2, h_bb);
    if (FeatComp::_feat_check("dR_l1_l2_boosted_htt_met")) feats["dR_l1_l2_boosted_htt_met"] = FeatComp::delta_r_boosted(l_1, l_2, h_tt_met);
    if (FeatComp::_feat_check("dR_l1_l2_boosted_sv"))      feats["dR_l1_l2_boosted_sv"]      = FeatComp::delta_r_boosted(l_1, l_2, svfit);
    
    // Masses
    if (FeatComp::_feat_check("h_tt_met_mt"))  feats["h_tt_met_mt"]  = FeatComp::calc_mt(h_tt_met, met);
    if (FeatComp::_feat_check("diH_mass_met")) feats["diH_mass_met"] = (h_bb+h_tt_met).M();
    if (FeatComp::_feat_check("diH_mass_sv"))  feats["diH_mass_sv"]  = (h_bb+svfit).M();
    if (FeatComp::_feat_check("diH_mass_vis")) feats["diH_mass_vis"] = (h_bb+h_tt_vis).M();
    if (FeatComp::_feat_check("diH_mass_X"))   feats["diH_mass_X"]   = (h_bb+h_tt_met).M()-h_tt_met.M()-h_bb.M();

    // Angle Deltas
    if (FeatComp::_feat_check("phi_met")) feats["phi_met"] = FeatComp::calc_phi(l_1, l_2, b_1, b_2, h_bb+h_tt_met);
    if (FeatComp::_feat_check("phi"))     feats["phi"]     = FeatComp::calc_phi(l_1, l_2, b_1, b_2, hh);
    if (FeatComp::_feat_check("phi1_met")) feats["phi1_met"] = FeatComp::calc_phi_1(l_1, l_2, h_tt_met, h_bb+h_tt_met);
    if (FeatComp::_feat_check("phi1"))     feats["phi1"]     = FeatComp::calc_phi_1(l_1, l_2, svfit, hh);
    if (FeatComp::_feat_check("phi2_met")) feats["phi2_met"] = FeatComp::calc_phi_1(b_1, b_2, h_tt_met, h_bb+h_tt_met);
    if (FeatComp::_feat_check("phi2"))     feats["phi2"]     = FeatComp::calc_phi_1(b_1, b_2, svfit, hh);
    if (FeatComp::_feat_check("costheta_star_met")) feats["costheta_star_met"] = FeatComp::calc_cos_delta_star(h_tt_met, h_bb+h_tt_met);
    if (FeatComp::_feat_check("costheta_star"))     feats["costheta_star"]     = FeatComp::calc_cos_delta_star(svfit, hh);
    if (FeatComp::_feat_check("costheta_l1_httmet"))      feats["costheta_l1_httmet"]      = FeatComp::calc_cos_delta(l_1, h_tt_met);
    if (FeatComp::_feat_check("costheta_l1_htt"))         feats["costheta_l1_httmet"]      = FeatComp::calc_cos_delta(l_1, svfit);
    if (FeatComp::_feat_check("costheta_l2_httmet"))      feats["costheta_l2_httmet"]      = FeatComp::calc_cos_delta(l_2, h_tt_met);
    if (FeatComp::_feat_check("costheta_l2_htt"))         feats["costheta_l2_httmet"]      = FeatComp::calc_cos_delta(l_2, svfit);
    if (FeatComp::_feat_check("costheta_met_httmet"))     feats["costheta_met_httmet"]     = FeatComp::calc_cos_delta(met, h_tt_met);
    if (FeatComp::_feat_check("costheta_met_htt"))        feats["costheta_met_htt"]        = FeatComp::calc_cos_delta(met, svfit);
    if (FeatComp::_feat_check("costheta_met_hbb"))        feats["costheta_met_hbb"]        = FeatComp::calc_cos_delta(met, h_bb);
    if (FeatComp::_feat_check("costheta_b1_hbb"))         feats["costheta_b1_hbb"]         = FeatComp::calc_cos_delta(b_1, h_bb);
    if (FeatComp::_feat_check("costheta_b2_hbb"))         feats["costheta_b2_hbb"]         = FeatComp::calc_cos_delta(b_2, h_bb);
    if (FeatComp::_feat_check("costheta_htt_hh_vis"))     feats["costheta_htt_hh_vis"]     = FeatComp::calc_cos_delta(svfit, h_tt_vis+h_bb);
    if (FeatComp::_feat_check("costheta_htt_met_hh_vis")) feats["costheta_htt_met_hh_vis"] = FeatComp::calc_cos_delta(h_tt_met, h_tt_vis+h_bb);
    if (FeatComp::_feat_check("costheta_hbb_hh_vis"))     feats["costheta_hbb_hh_vis"]     = FeatComp::calc_cos_delta(h_bb, h_tt_vis+h_bb);
    if (FeatComp::_feat_check("costheta_htt_hh_met"))     feats["costheta_htt_hh_met"]     = FeatComp::calc_cos_delta(svfit, h_tt_met+h_bb);
    if (FeatComp::_feat_check("costheta_htt_met_hh_met")) feats["costheta_htt_met_hh_met"] = FeatComp::calc_cos_delta(h_tt_met, h_tt_met+h_bb);
    if (FeatComp::_feat_check("costheta_hbb_hh_met"))     feats["costheta_hbb_hh_met"]     = FeatComp::calc_cos_delta(h_bb, h_tt_met+h_bb);
    if (FeatComp::_feat_check("costheta_htt_hh"))         feats["costheta_htt_hh_vis"]     = FeatComp::calc_cos_delta(svfit, hh);
    if (FeatComp::_feat_check("costheta_htt_met_hh"))     feats["costheta_htt_met_hh_vis"] = FeatComp::calc_cos_delta(h_tt_met, hh);
    if (FeatComp::_feat_check("costheta_hbb_hh"))         feats["costheta_hbb_hh_vis"]     = FeatComp::calc_cos_delta(h_bb, hh);            

    return feats;
}

inline bool FeatComp::_feat_check(std::string feat) {return (_all ? true : std::find(_requested.begin(), _requested.end(), feat) != _requested.end())}

inline float FeatComp::delta_phi(const TLorentzVector& v_0, const TLorentzVector& v_1) return std::abs(v_0.DeltaPhi(v_1));

inline float FeatComp::delta_eta(const TLorentzVector& v_0, const TLorentzVector& v_1) return std::abs(v_0.Eta()-v_1.Eta());

inline float  FeatComp::delta_r_boosted(const TLorentzVector& v_0, const TLorentzVector& v_1, const TLorentzVector& ref){
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    using namespace ROOT::Math::VectorUtil;
    return DeltaR(boost(v_0, ref.BoostToCM()), boost(v_1, ref.BoostToCM()));
}

inline float FeatComp::calc_mt(const TLorentzVector& v, const TLorentzVector& met) return std::sqrt(2.0*v.Pt()*met.Pt()*(1.0-std::cos(v.DeltaPhi(met))));

inline float FeatComp::calc_phi(const TLorentzVector& l_1, const TLorentzVector& l_2,
                                const TLorentzVector& b_1, const TLorentzVector& b_2, const TLorentzVector& hh) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    using namespace ROOT::Math::VectorUtil;
    return Angle(boost(l_1, hh.BoostToCM()).Vect().Cross(boost(l_2, hh.BoostToCM()).Vect()),
                 boost(b_1, hh.BoostToCM()).Vect().Cross(boost(b_2, hh.BoostToCM()).Vect()));
}

inline float FeatComp::calc_phi_1(const TLorentzVector& v_0, const TLorentzVector& v_1, const TLorentzVector& h, const TLorentzVector& hh) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    
    using namespace ROOT::Math::VectorUtil;
    return Angle(boost(v_0, hh.BoostToCM()).Vect().Cross(boost(v_1, hh.BoostToCM()).Vect()),
                 boost(h, hh.BoostToCM()).Vect().Cross(ROOT::Math::Cartesian3D<>(0, 0, 1)));
}

inline float FeatComp::calc_cos_delta_star(const TLorentzVector& v, const TLorentzVector& hh) { 
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    
    using namespace ROOT::Math::VectorUtil;
    return CosTheta(boost(v, hh.BoostToCM()), ROOT::Math::Cartesian3D<>(0, 0, 1));
}

inline float FeatComp::calc_cos_delta(const TLorentzVector& v, const TLorentzVector& r) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    
    using namespace ROOT::Math::VectorUtil;
    return  CosTheta(boost(v, r.BoostToCM()), r);
}

void FeatComp::_add_jet_flags(const float& b_1_csv, const float& b_2_csv,
                              const float& b_1_deepcsv, const float& b_2_deepcsv,
                              std::map<std::string, float>& feats) {
    int tag_1(0), tag_2(0);
    float csv_1(_use_deep_csv ? b_1_deepcsv : b_1_csv), csv_2(_use_deep_csv ? b_2_deepcsv : b_2_csv);
    for (float wp : (_use_deep_csv ? _deep_csv_wps : _csv_wps)) {
        if (csv_1 >= wp) tag_1++;
        if (csv_2 >= wp) tag_2++;
    }
    feats["jet_1_quality"] = tag_1;
    feats["jet_2_quality"] = tag_2;
}
