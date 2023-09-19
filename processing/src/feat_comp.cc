#include "cms_hh_proc_interface/processing/interface/feat_comp.hh"

FeatComp::FeatComp(bool return_all, std::vector<std::string> requested, bool use_deep_bjet_wps) {
    _all = return_all;
    _requested = requested;
    _use_deep_bjet_wps = use_deep_bjet_wps;
}

FeatComp::~FeatComp() {}

std::map<std::string, float> FeatComp::process(const LorentzVector& b_1,
										       const LorentzVector& b_2, 
										       const LorentzVector& l_1,
										       const LorentzVector& l_2,
										       const LorentzVector& met,
										       const LorentzVector& svfit,
                                               const LorentzVector& vbf_1,
                                               const LorentzVector& vbf_2,
                                               const LorentzVector& Nu_1,
                                               const LorentzVector& Nu_2,
										       const float& hh_kinfit_m,
                                               const bool& is_boosted,
                                               const float& b_1_csv,
                                               const float& b_2_csv,
                                               Channel channel,
                                               Year year,
                                               const int& n_vbf,
                                               const float& dau1_flav,
                                               const float& dau2_flav,
                                               const bool& svfit_conv,
                                               const bool& hh_kinfit_conv,
                                               const float& b_1_cvsl,
										       const float& b_2_cvsl,
										       const float& vbf_1_cvsl,
										       const float& vbf_2_cvsl,
										       const float& b_1_cvsb,
										       const float& b_2_cvsb,
										       const float& vbf_1_cvsb,
										       const float& vbf_2_cvsb,
                                               const float& met_et,
                                               const float& met_phi,
                                               const float& DeepMET_ResponseTune_px,
                                               const float& DeepMET_ResponseTune_py,
                                               const float& DeepMET_ResolutionTune_px,
                                               const float& DeepMET_ResolutionTune_py) {
    /* Compute HL features from base event*/

    bool use_vbf = n_vbf >= 2;
    // Extra vectors
    LorentzVector h_bb(b_1 + b_2);
    LorentzVector h_tt_vis(l_1 + l_2);
    LorentzVector h_tt_met(h_tt_vis + met);
    LorentzVector hh(h_bb.Px()+svfit.Px(), h_bb.Py()+svfit.Py(), h_bb.Pz()+svfit.Pz(), hh_kinfit_m);  // I assume 4th component is a mass, but who knows...
    if (!hh_kinfit_conv) {  // HHKinFit didn't converge
        hh = svfit_conv ? h_bb+svfit : h_bb+h_tt_met;
    } else if (!svfit_conv) {  // HHKinFit converge but SVFit didn't
        hh = LorentzVector(h_bb.Px()+h_tt_met.Px(), h_bb.Py()+h_tt_met.Py(), h_bb.Pz()+h_tt_met.Pz(), hh_kinfit_m);
    }
    std::map<std::string, float> feats;

    // Categoricals
    if (FeatComp::_feat_check("boosted"))    feats["boosted"]    = is_boosted;
    if (FeatComp::_feat_check("is_vbf"))     feats["is_vbf"]     = n_vbf >= 2;
    if (FeatComp::_feat_check("channel"))    feats["channel"]    = channel;
    if (FeatComp::_feat_check("year"))       feats["year"]       = year;
    if (FeatComp::_feat_check("svfit_conv")) feats["svfit_conv"] = svfit_conv;
    if (FeatComp::_feat_check("n_vbf"))      feats["n_vbf"]      = n_vbf;
    FeatComp::_add_btag_flags(year, b_1_csv, b_2_csv, feats);

    // continuous 

    // btag csvl and csvb scores
    if (FeatComp::_feat_check("b_1_csv"))     feats["b_1_csv"]    = b_1_cvsl;
    if (FeatComp::_feat_check("b_2_csv"))     feats["b_2_csv"]    = b_2_cvsl;
    if (FeatComp::_feat_check("b_1_cvsl"))    feats["b_1_cvsl"]   = b_1_cvsl;
    if (FeatComp::_feat_check("b_2_cvsl"))    feats["b_2_cvsl"]   = b_2_cvsl;
    if (FeatComp::_feat_check("vbf_1_cvsl"))  feats["vbf_1_cvsl"] = use_vbf ? vbf_1_cvsl : -1;
    if (FeatComp::_feat_check("vbf_2_cvsl"))  feats["vbf_2_cvsl"] = use_vbf ? vbf_2_cvsl : -1;
    if (FeatComp::_feat_check("b_1_cvsb"))    feats["b_1_cvsb"]   = b_1_cvsb;
    if (FeatComp::_feat_check("b_2_cvsb"))    feats["b_2_cvsb"]   = b_2_cvsb;
    if (FeatComp::_feat_check("vbf_1_cvsb"))  feats["vbf_1_cvsb"] = use_vbf ? vbf_1_cvsb : -1;
    if (FeatComp::_feat_check("vbf_2_cvsb"))  feats["vbf_2_cvsb"] = use_vbf ? vbf_2_cvsb : -1;
    // Delta phi
    if (FeatComp::_feat_check("dphi_l1_l2"))      feats["dphi_l1_l2"]      = FeatComp::delta_phi(l_1, l_2);
    if (FeatComp::_feat_check("dphi_b1_b2"))      feats["dphi_b1_b2"]      = FeatComp::delta_phi(b_1, b_2);
    if (FeatComp::_feat_check("dphi_l1_met"))     feats["dphi_l1_met"]     = FeatComp::delta_phi(l_1, met);
    if (FeatComp::_feat_check("dphi_l2_met"))     feats["dphi_l2_met"]     = FeatComp::delta_phi(l_2, met);
    if (FeatComp::_feat_check("dphi_l1l2_met"))   feats["dphi_l1l2_met"]   = FeatComp::delta_phi(h_tt_vis, met);
    if (FeatComp::_feat_check("dphi_sv_met"))     feats["dphi_sv_met"]     = svfit_conv ? FeatComp::delta_phi(svfit, met) : -1;
    if (FeatComp::_feat_check("dphi_httmet_met")) feats["dphi_httmet_met"] = FeatComp::delta_phi(h_tt_met, met);
    if (FeatComp::_feat_check("dphi_httvis_met")) feats["dphi_httvis_met"] = FeatComp::delta_phi(h_tt_vis, met);
    if (FeatComp::_feat_check("dphi_hbb_met"))    feats["dphi_hbb_met"]    = FeatComp::delta_phi(h_bb, met);
    if (FeatComp::_feat_check("dphi_hbb_sv"))     feats["dphi_hbb_sv"]     = svfit_conv ? FeatComp::delta_phi(h_bb, svfit) :-1;
    if (FeatComp::_feat_check("dphi_hbb_httmet")) feats["dphi_hbb_httmet"] = FeatComp::delta_phi(h_bb, h_tt_met);
    if (FeatComp::_feat_check("dphi_vbf1_vbf2"))  feats["dphi_vbf1_vbf2"]  = use_vbf ? FeatComp::delta_phi(vbf_1, vbf_2) :-1;
    if (FeatComp::_feat_check("dphi_vbf1_met"))   feats["dphi_vbf1_met"]   = use_vbf ? FeatComp::delta_phi(vbf_1, met)   :-1;
    if (FeatComp::_feat_check("dphi_vbf2_met"))   feats["dphi_vbf2_met"]   = use_vbf ? FeatComp::delta_phi(vbf_2, met)   :-1;

    // Delta eta
    if (FeatComp::_feat_check("deta_l1_l2"))      feats["deta_l1_l2"]      = FeatComp::delta_eta(l_1, l_2);
    if (FeatComp::_feat_check("deta_b1_b2"))      feats["deta_b1_b2"]      = FeatComp::delta_eta(b_1, b_2);
    if (FeatComp::_feat_check("deta_hbb_sv"))     feats["deta_hbb_sv"]     = svfit_conv ? FeatComp::delta_eta(h_bb, svfit) : -1;
    if (FeatComp::_feat_check("deta_hbb_httmet")) feats["deta_hbb_httmet"] = FeatComp::delta_eta(h_bb, h_tt_met);
    if (FeatComp::_feat_check("deta_hbb_httvis")) feats["deta_hbb_httvis"] = FeatComp::delta_eta(h_bb, h_tt_vis);
    if (FeatComp::_feat_check("deta_vbf1_vbf2"))  feats["deta_vbf1_vbf2"]  = use_vbf ? FeatComp::delta_eta(vbf_1, vbf_2) : -1;

    // Delta R
    if (FeatComp::_feat_check("dR_l1_l2"))      feats["dR_l1_l2"]      = FeatComp::delta_r(l_1, l_2);
    if (FeatComp::_feat_check("dR_b1_b2"))      feats["dR_b1_b2"]      = FeatComp::delta_r(b_1, b_2);
    if (FeatComp::_feat_check("dR_hbb_httmet")) feats["dR_hbb_httmet"] = FeatComp::delta_r(h_bb, h_tt_met);
    if (FeatComp::_feat_check("dR_hbb_sv"))     feats["dR_hbb_sv"]     = svfit_conv ? FeatComp::delta_r(h_bb, svfit) :-1;
    if (FeatComp::_feat_check("dR_vbf1_vbf2"))  feats["dR_vbf1_vbf2"]  = use_vbf ? FeatComp::delta_r(vbf_1, vbf_2) :-1;

    if (FeatComp::_feat_check("dR_b1_b2_x_h_bb_pT"))     feats["dR_b1_b2_x_h_bb_pT"]     = FeatComp::delta_r(b_1, b_2)*h_bb.Pt();
    if (FeatComp::_feat_check("dR_l1_l2_x_h_tt_met_pT")) feats["dR_l1_l2_x_h_tt_met_pT"] = FeatComp::delta_r(l_1, l_2)*h_tt_met.Pt();
    if (FeatComp::_feat_check("dR_l1_l2_x_sv_pT"))       feats["dR_l1_l2_x_sv_pT"]       = svfit_conv ? FeatComp::delta_r(l_1, l_2)*svfit.Pt() :-1;
    if (FeatComp::_feat_check("dR_b1_b2_boosted_hbb"))     feats["dR_b1_b2_boosted_hbb"]     = FeatComp::delta_r_boosted(b_1, b_2, h_bb);
    if (FeatComp::_feat_check("dR_l1_l2_boosted_htt_met")) feats["dR_l1_l2_boosted_htt_met"] = FeatComp::delta_r_boosted(l_1, l_2, h_tt_met);
    if (FeatComp::_feat_check("dR_l1_l2_boosted_sv"))      feats["dR_l1_l2_boosted_sv"]      = svfit_conv ? FeatComp::delta_r_boosted(l_1, l_2, svfit) :-1;
    if (FeatComp::_feat_check("min_dR_vbfj_l")) 
    {
      float min_dR_vbfj_l = FeatComp::delta_r(vbf_1, l_1); 
      float tmp_dr = FeatComp::delta_r(vbf_1, l_2);
      if (tmp_dr < min_dR_vbfj_l) min_dR_vbfj_l = tmp_dr;
      tmp_dr = FeatComp::delta_r(vbf_2, l_1);
      if (tmp_dr < min_dR_vbfj_l) min_dR_vbfj_l = tmp_dr;
      tmp_dr = FeatComp::delta_r(vbf_2, l_2);
      if (tmp_dr < min_dR_vbfj_l) min_dR_vbfj_l = tmp_dr;
      feats["min_dR_vbfj_b"]      = min_dR_vbfj_l;
    }
    if (FeatComp::_feat_check("min_dR_vbfj_b")) 
    {
      float min_dR_vbfj_b = FeatComp::delta_r(vbf_1, b_1); 
      float tmp_dr = FeatComp::delta_r(vbf_1, b_2);
      if (tmp_dr < min_dR_vbfj_b) min_dR_vbfj_b = tmp_dr;
      tmp_dr = FeatComp::delta_r(vbf_2, b_1);
      if (tmp_dr < min_dR_vbfj_b) min_dR_vbfj_b = tmp_dr;
      tmp_dr = FeatComp::delta_r(vbf_2, b_2);
      if (tmp_dr < min_dR_vbfj_b) min_dR_vbfj_b = tmp_dr;
      feats["min_dR_vbfj_b"]      = min_dR_vbfj_b;
    }

    // Masses
    if (FeatComp::_feat_check("sv_mass"))       feats["sv_mass"]       = svfit_conv ? svfit.M() :-1;
    if (FeatComp::_feat_check("h_tt_vis_mass")) feats["h_tt_vis_mass"] = h_tt_vis.M();
    if (FeatComp::_feat_check("h_bb_mass"))     feats["h_bb_mass"]     = h_bb.M();
    if (FeatComp::_feat_check("hh_kinfit_m"))   feats["hh_kinfit_m"]   = hh_kinfit_conv ? hh_kinfit_m : 0;
    if (FeatComp::_feat_check("sv_mt"))         feats["sv_mt"]         = svfit_conv ? FeatComp::calc_mt(svfit, met) :-1;
    if (FeatComp::_feat_check("h_tt_met_mt"))   feats["h_tt_met_mt"]   = FeatComp::calc_mt(h_tt_met, met);
    if (FeatComp::_feat_check("l_1_mt"))        feats["l_1_mt"]        = FeatComp::calc_mt(l_1, met);
    if (FeatComp::_feat_check("l_2_mt"))        feats["l_2_mt"]        = FeatComp::calc_mt(l_2, met);
    if (FeatComp::_feat_check("ll_mt"))         feats["ll_mt"]         = FeatComp::calc_mt(l_1, l_2);
    if (FeatComp::_feat_check("mt_tot"))        feats["mt_tot"]        = FeatComp::calc_mt_tot(l_1, l_2, met);
    if (FeatComp::_feat_check("diH_mass_met"))  feats["diH_mass_met"]  = (h_bb+h_tt_met).M();
    if (FeatComp::_feat_check("diH_mass_sv"))   feats["diH_mass_sv"]   = svfit_conv ? (h_bb+svfit).M() :-1;
    if (FeatComp::_feat_check("diH_mass_vis"))  feats["diH_mass_vis"]  = (h_bb+h_tt_vis).M();
    if (FeatComp::_feat_check("diH_mass_X"))    feats["diH_mass_X"]    = (h_bb+h_tt_met).M()-h_tt_met.M()-h_bb.M();
    if (FeatComp::_feat_check("vbf_invmass"))   feats["vbf_invmass"]   = use_vbf ? (vbf_1+vbf_2).M() :-1;

    // Momenta
    if (FeatComp::_feat_check("sv_pT"))   feats["sv_pT"]   = svfit_conv ? svfit.Pt() : -1;
    if (FeatComp::_feat_check("h_bb_pT")) feats["h_bb_pT"] = h_bb.Pt();
    if (FeatComp::_feat_check("hh_pT"))   feats["hh_pT"]   = hh.Pt();

    // Energies
    if (FeatComp::_feat_check("sv_E"))   feats["sv_E"]   = svfit_conv ? svfit.E() :-1;
    if (FeatComp::_feat_check("h_bb_E")) feats["h_bb_E"] = h_bb.E();
    if (FeatComp::_feat_check("hh_E"))   feats["hh_E"]   = hh.E();

    // Angle Deltas
    if (FeatComp::_feat_check("phi_met")) feats["phi_met"] = FeatComp::calc_phi(l_1, l_2, b_1, b_2, h_bb+h_tt_met);
    if (FeatComp::_feat_check("phi"))     feats["phi"]     = FeatComp::calc_phi(l_1, l_2, b_1, b_2, hh);
    if (FeatComp::_feat_check("phi1_met")) feats["phi1_met"] = FeatComp::calc_phi_1(l_1, l_2, h_tt_met, h_bb+h_tt_met);
    if (FeatComp::_feat_check("phi1"))     feats["phi1"]     = svfit_conv ? FeatComp::calc_phi_1(l_1, l_2, svfit, hh) : -1;
    if (FeatComp::_feat_check("phi2_met")) feats["phi2_met"] = FeatComp::calc_phi_1(b_1, b_2, h_tt_met, h_bb+h_tt_met);
    if (FeatComp::_feat_check("phi2"))     feats["phi2"]     = svfit_conv ? FeatComp::calc_phi_1(b_1, b_2, svfit, hh) : -1;
    if (FeatComp::_feat_check("costheta_star_met")) feats["costheta_star_met"] = FeatComp::calc_cos_delta_star(h_tt_met, h_bb+h_tt_met);
    if (FeatComp::_feat_check("costheta_star"))     feats["costheta_star"]     = svfit_conv ? FeatComp::calc_cos_delta_star(svfit, hh) : -2;
    if (FeatComp::_feat_check("costheta_l1_httmet"))      feats["costheta_l1_httmet"]      = FeatComp::calc_cos_delta(l_1, h_tt_met);
    if (FeatComp::_feat_check("costheta_l1_htt"))         feats["costheta_l1_htt"]         = svfit_conv ? FeatComp::calc_cos_delta(l_1, svfit) :-2;
    if (FeatComp::_feat_check("costheta_l2_httmet"))      feats["costheta_l2_httmet"]      = FeatComp::calc_cos_delta(l_2, h_tt_met);
    if (FeatComp::_feat_check("costheta_l2_htt"))         feats["costheta_l2_htt"]         = svfit_conv ? FeatComp::calc_cos_delta(l_2, svfit) :-2;
    if (FeatComp::_feat_check("costheta_met_httmet"))     feats["costheta_met_httmet"]     = FeatComp::calc_cos_delta(met, h_tt_met);
    if (FeatComp::_feat_check("costheta_met_htt"))        feats["costheta_met_htt"]        = svfit_conv ? FeatComp::calc_cos_delta(met, svfit) :-2;
    if (FeatComp::_feat_check("costheta_met_hbb"))        feats["costheta_met_hbb"]        = FeatComp::calc_cos_delta(met, h_bb);
    if (FeatComp::_feat_check("costheta_b1_hbb"))         feats["costheta_b1_hbb"]         = FeatComp::calc_cos_delta(b_1, h_bb);
    if (FeatComp::_feat_check("costheta_htt_hh_vis"))     feats["costheta_htt_hh_vis"]     = svfit_conv ? FeatComp::calc_cos_delta(svfit, h_tt_vis+h_bb) :-2;
    if (FeatComp::_feat_check("costheta_htt_met_hh_vis")) feats["costheta_htt_met_hh_vis"] = FeatComp::calc_cos_delta(h_tt_met, h_tt_vis+h_bb);
    if (FeatComp::_feat_check("costheta_hbb_hh_vis"))     feats["costheta_hbb_hh_vis"]     = FeatComp::calc_cos_delta(h_bb, h_tt_vis+h_bb);
    if (FeatComp::_feat_check("costheta_htt_hh_met"))     feats["costheta_htt_hh_met"]     = svfit_conv ? FeatComp::calc_cos_delta(svfit, h_tt_met+h_bb) :-2;
    if (FeatComp::_feat_check("costheta_htt_met_hh_met")) feats["costheta_htt_met_hh_met"] = FeatComp::calc_cos_delta(h_tt_met, h_tt_met+h_bb);
    if (FeatComp::_feat_check("costheta_hbb_hh_met"))     feats["costheta_hbb_hh_met"]     = FeatComp::calc_cos_delta(h_bb, h_tt_met+h_bb);
    if (FeatComp::_feat_check("costheta_htt_hh"))         feats["costheta_htt_hh"]         = svfit_conv ? FeatComp::calc_cos_delta(svfit, hh) :-2;
    if (FeatComp::_feat_check("costheta_htt_met_hh"))     feats["costheta_htt_met_hh"]     = FeatComp::calc_cos_delta(h_tt_met, hh);
    if (FeatComp::_feat_check("costheta_hbb_hh"))         feats["costheta_hbb_hh"]         = FeatComp::calc_cos_delta(h_bb, hh);

    //Centralities
    if (FeatComp::_feat_check("l_1_centrality"))      feats["l_1_centrality"]      =  use_vbf ? FeatComp::calc_centrality(l_1, vbf_1, vbf_2) :-1;
    if (FeatComp::_feat_check("l_2_centrality"))      feats["l_2_centrality"]      =  use_vbf ? FeatComp::calc_centrality(l_2, vbf_1, vbf_2) :-1;
    if (FeatComp::_feat_check("b_1_centrality"))      feats["b_1_centrality"]      =  use_vbf ? FeatComp::calc_centrality(b_1, vbf_1, vbf_2) :-1;
    if (FeatComp::_feat_check("b_2_centrality"))      feats["b_2_centrality"]      =  use_vbf ? FeatComp::calc_centrality(b_2, vbf_1, vbf_2) :-1;
    if (FeatComp::_feat_check("h_bb_centrality"))     feats["h_bb_centrality"]     =  use_vbf ? FeatComp::calc_centrality(h_bb, vbf_1, vbf_2) :-1;
    if (FeatComp::_feat_check("h_tt_vis_centrality")) feats["h_tt_vis_centrality"] =  use_vbf ? FeatComp::calc_centrality(h_tt_vis, vbf_1, vbf_2) :-1;
    if (FeatComp::_feat_check("hh_centrality"))       feats["hh_centrality"]       =  use_vbf ? FeatComp::calc_hh_centrality(h_bb, h_tt_vis, vbf_1, vbf_2) :-1;

    // Assorted VBF
    if (FeatComp::_feat_check("vbf_eta_prod")) feats["vbf_eta_prod"] = use_vbf ? vbf_1.eta()*vbf_2.eta() : std::nanf("1");

    // Htautau inputs
    float DeepMET_ResolutionTune_phi = std::atan2(DeepMET_ResponseTune_py, DeepMET_ResponseTune_px);
    if (FeatComp::_feat_check("dmet_resp_px"))      feats["dmet_resp_px"] = FeatComp::calc_dmet_px(DeepMET_ResponseTune_px, DeepMET_ResponseTune_py, DeepMET_ResolutionTune_phi);
    if (FeatComp::_feat_check("dmet_resp_py"))      feats["dmet_resp_py"] = FeatComp::calc_dmet_py(DeepMET_ResponseTune_px, DeepMET_ResponseTune_py, DeepMET_ResolutionTune_phi);
    if (FeatComp::_feat_check("dmet_reso_px"))      feats["dmet_reso_px"] = FeatComp::calc_dmet_px(DeepMET_ResolutionTune_px,DeepMET_ResolutionTune_py,DeepMET_ResolutionTune_phi);
    if (FeatComp::_feat_check("dmet_reso_py"))      feats["dmet_reso_py"] = FeatComp::calc_dmet_py(DeepMET_ResolutionTune_px, DeepMET_ResolutionTune_py, DeepMET_ResolutionTune_phi);
    if (FeatComp::_feat_check("met_px"))            feats["met_px"] = met_et * std::cos(met_phi);
    if (FeatComp::_feat_check("met_py"))            feats["met_py"] = met_et * std::sin(met_phi);
    //if (FeatComp::_feat_check("dau1_dphi"))         feats["dau1_dphi"] = FeatComp::mpi_to_pi(l_1.Phi() - DeepMET_ResolutionTune_phi);
    //if (FeatComp::_feat_check("dau2_dphi"))         feats["dau2_dphi"] = FeatComp::mpi_to_pi(l_2.Phi() - DeepMET_ResolutionTune_phi);
    if (FeatComp::_feat_check("dau1_charge"))       feats["dau1_charge"] = dau1_flav/std::abs(dau1_flav);
    if (FeatComp::_feat_check("dau2_charge"))       feats["dau2_charge"] = dau2_flav/std::abs(dau2_flav);
    //if (FeatComp::_feat_check("genNu1_dphi"))       feats["genNu1_dphi"] = FeatComp::mpi_to_pi(Nu_1.Phi() - DeepMET_ResolutionTune_phi);
    //if (FeatComp::_feat_check("genNu2_dphi"))       feats["genNu2_dphi"] = FeatComp::mpi_to_pi(Nu_2.Phi() - DeepMET_ResolutionTune_phi);
    if (FeatComp::_feat_check("dau1_px"))           feats["dau1_px"] = l_1.Px();
    if (FeatComp::_feat_check("dau1_py"))           feats["dau1_py"] = l_1.Py();
    if (FeatComp::_feat_check("dau1_pz"))           feats["dau1_pz"] = l_1.Pz();
    if (FeatComp::_feat_check("dau1_m"))            feats["dau1_m"] = l_1.M();
    if (FeatComp::_feat_check("dau2_px"))           feats["dau2_px"] = l_2.Px();
    if (FeatComp::_feat_check("dau2_py"))           feats["dau2_py"] = l_2.Py();
    if (FeatComp::_feat_check("dau2_pz"))           feats["dau2_pz"] = l_2.Pz();
    if (FeatComp::_feat_check("dau2_m"))            feats["dau2_m"] = l_2.M();
    //if (FeatComp::_feat_check("ditau_deltaphi"))    feats["ditau_deltaphi"] = std::abs(FeatComp::mpi_to_pi(FeatComp::mpi_to_pi(l_1.Phi() - DeepMET_ResolutionTune_phi) - FeatComp::mpi_to_pi(l_2.Phi() - DeepMET_ResolutionTune_phi)));
    //if (FeatComp::_feat_check("ditau_deltaeta"))    feats["ditau_deltaeta"] = std::abs(l_1.Eta() - l_2.Eta());
    if (FeatComp::_feat_check("genNu1_px"))         feats["genNu1_px"] = Nu_1.Px();
    if (FeatComp::_feat_check("genNu1_py"))         feats["genNu1_py"] = Nu_1.Py();
    if (FeatComp::_feat_check("genNu1_pz"))         feats["genNu1_pz"] = Nu_1.Pz();
    if (FeatComp::_feat_check("genNu2_px"))         feats["genNu2_px"] = Nu_2.Px();
    if (FeatComp::_feat_check("genNu2_py"))         feats["genNu2_py"] = Nu_2.Py();
    if (FeatComp::_feat_check("genNu2_pz"))         feats["genNu2_pz"] = Nu_2.Pz();
    //if (FeatComp::_feat_check("bjet1_dphi"))        feats["bjet1_dphi"] = std::abs(FeatComp::mpi_to_pi(b_1.Phi() - DeepMET_ResolutionTune_phi));
    if (FeatComp::_feat_check("bjet1_px"))          feats["bjet1_px"] = b_1.Px();
    if (FeatComp::_feat_check("bjet1_py"))          feats["bjet1_py"] = b_1.Py();
    if (FeatComp::_feat_check("bjet1_pz"))          feats["bjet1_pz"] = b_1.Pz();
    //if (FeatComp::_feat_check("bjet2_dphi"))        feats["bjet2_dphi"] = std::abs(FeatComp::mpi_to_pi(b_2.Phi() - DeepMET_ResolutionTune_phi));
    if (FeatComp::_feat_check("bjet2_px"))          feats["bjet2_px"] = b_2.Px();
    if (FeatComp::_feat_check("bjet2_py"))          feats["bjet2_py"] = b_2.Py();
    if (FeatComp::_feat_check("bjet2_pz"))          feats["bjet2_pz"] = b_2.Pz();


    // Assorted HL 
    if (FeatComp::_feat_check("p_zetavisible")) feats["p_zetavisible"] = FeatComp::calc_pzeta_visible(l_1, l_2);
    if (FeatComp::_feat_check("p_zeta"))        feats["p_zeta"]        = FeatComp::calc_pzeta(l_1, l_2, met);
    if (FeatComp::_feat_check("top_1_mass") || FeatComp::_feat_check("top_2_mass")) {
        std::pair<float, float> top_masses = FeatComp::calc_top_masses(l_1, l_2, b_1, b_2, met);
        if (FeatComp::_feat_check("top_1_mass")) feats["top_1_mass"] = top_masses.first;
        if (FeatComp::_feat_check("top_2_mass")) feats["top_2_mass"] = top_masses.second;
    }
    return feats;
}

inline bool FeatComp::_feat_check(std::string feat) {return (_all ? true : std::find(_requested.begin(), _requested.end(), feat) != _requested.end());}

inline float FeatComp::delta_phi(const LorentzVector& v_0, const LorentzVector& v_1) {return std::abs(ROOT::Math::VectorUtil::DeltaPhi(v_0, v_1));}

inline float FeatComp::delta_eta(const LorentzVector& v_0, const LorentzVector& v_1) {return std::abs(v_0.Eta()-v_1.Eta());}

inline float FeatComp::delta_r(const LorentzVector& v_0, const LorentzVector& v_1) {return ROOT::Math::VectorUtil::DeltaR(v_0, v_1);}

inline float  FeatComp::delta_r_boosted(const LorentzVector& v_0, const LorentzVector& v_1, const LorentzVector& ref){
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    using namespace ROOT::Math::VectorUtil;
    return DeltaR(boost(v_0, ref.BoostToCM()), boost(v_1, ref.BoostToCM()));
}


inline float FeatComp::calc_mt(const LorentzVector& v, const LorentzVector& met) {
    return std::sqrt(2.0*v.Pt()*met.Pt()*(1.0-std::cos(ROOT::Math::VectorUtil::DeltaPhi(v,met))));
}

inline float FeatComp::calc_phi(const LorentzVector& l_1, const LorentzVector& l_2,
                                const LorentzVector& b_1, const LorentzVector& b_2, const LorentzVector& hh) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    using namespace ROOT::Math::VectorUtil;
    return Angle(boost(l_1, hh.BoostToCM()).Vect().Cross(boost(l_2, hh.BoostToCM()).Vect()),
                 boost(b_1, hh.BoostToCM()).Vect().Cross(boost(b_2, hh.BoostToCM()).Vect()));
}

inline float FeatComp::calc_phi_1(const LorentzVector& v_0, const LorentzVector& v_1, const LorentzVector& h, const LorentzVector& hh) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    
    using namespace ROOT::Math::VectorUtil;
    return Angle(boost(v_0, hh.BoostToCM()).Vect().Cross(boost(v_1, hh.BoostToCM()).Vect()),
                 boost(h, hh.BoostToCM()).Vect().Cross(ROOT::Math::Cartesian3D<>(0, 0, 1)));
}

inline float FeatComp::calc_cos_delta_star(const LorentzVector& v, const LorentzVector& hh) { 
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    
    using namespace ROOT::Math::VectorUtil;
    return CosTheta(boost(v, hh.BoostToCM()), ROOT::Math::Cartesian3D<>(0, 0, 1));
}

inline float FeatComp::calc_cos_delta(const LorentzVector& v, const LorentzVector& r) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */
    
    using namespace ROOT::Math::VectorUtil;
    return  CosTheta(boost(v, r.BoostToCM()), r);
}

inline float FeatComp::calc_dmet_px(const float& px,const float& py,const float& phi) {return std::cos(-phi)*px - std::sin(-phi)*py;}
inline float FeatComp::calc_dmet_py(const float& px,const float& py,const float& phi) {return std::sin(-phi)*px + std::cos(-phi)*py;}
float FeatComp::mpi_to_pi(const float& phi)
{
    float mod_phi(0);
    if (phi>TMath::Pi())
    {
        mod_phi = phi - 2*TMath::Pi();
        return mod_phi;
    }
    else if (phi<TMath::Pi())
    {
        mod_phi = phi + 2*TMath::Pi();
        return mod_phi;
    }
    else
    {
        return phi;
    }
}

void FeatComp::_add_btag_flags(Year year, const float& b_1_csv, const float& b_2_csv, std::map<std::string, float>& feats) {
    int tag_1(0), tag_2(0);
    for (float wp : (_use_deep_bjet_wps ? _deep_bjet_wps : _bjet_wps)[year]) {
        if (b_1_csv >= wp) tag_1++;
        if (b_2_csv >= wp) tag_2++;
    }
    if (FeatComp::_feat_check("jet_1_quality")) feats["jet_1_quality"] = tag_1;
    if (FeatComp::_feat_check("jet_2_quality")) feats["jet_2_quality"] = tag_2;
}

int FeatComp::_get_cvsl_flag(Year year, const float& score) {
    int tag(0);
    for (float wp : _cvsl_wps[year]) if (score >= wp) tag++;
    return tag;
}

int FeatComp::_get_cvsb_flag(Year year, const float& score) {
    int tag(0);
    for (float wp : _cvsb_wps[year]) if (score <= wp) tag++;
    return tag;
}

inline float FeatComp::calc_centrality( const LorentzVector& v, const LorentzVector& VBFjet1, const LorentzVector& VBFjet2)
{
    return std::abs(v.Eta() - 0.5*(VBFjet1.Eta() + VBFjet2.Eta() ))/(std::fabs(VBFjet1.Eta() - VBFjet2.Eta()));
}

inline float FeatComp::calc_hh_centrality( const LorentzVector& bh, const LorentzVector& tauh, const LorentzVector& VBFjet1, const LorentzVector& VBFjet2)
{
    float tmp_eta_minus = FeatComp::calcDeltaEtaMinus(bh, tauh, VBFjet1, VBFjet2);
    float tmp_eta_plus = FeatComp::calcDeltaEtaPlus(bh, tauh, VBFjet1, VBFjet2);
    if (tmp_eta_minus < tmp_eta_plus) return tmp_eta_minus;
    else return tmp_eta_plus;  
}

inline float FeatComp::calcDeltaEtaMinus(const LorentzVector& bh, const LorentzVector& tauh, const LorentzVector& VBFjet1, const LorentzVector& VBFjet2)  
{
    float bh_eta = bh.Eta(), tauh_eta = tauh.Eta(), VBFjet1_eta=VBFjet1.Eta(), VBFjet2_eta=VBFjet2.Eta(), min_eta_h, min_eta_vbf;
    if (bh_eta < tauh_eta) min_eta_h = bh_eta; 
    else min_eta_h = tauh_eta;
    if (VBFjet1_eta < VBFjet2_eta) min_eta_vbf = VBFjet1_eta; 
    else min_eta_vbf = VBFjet2_eta;
    return min_eta_h - min_eta_vbf;
}

inline float FeatComp::calcDeltaEtaPlus(const LorentzVector& bh, const LorentzVector& tauh, const LorentzVector& VBFjet1, const LorentzVector& VBFjet2)  
{
    float bh_eta = bh.Eta(), tauh_eta = tauh.Eta(), VBFjet1_eta=VBFjet1.Eta(), VBFjet2_eta=VBFjet2.Eta(), max_eta_h, max_eta_vbf;
    if (bh_eta > tauh_eta) max_eta_h = bh_eta; 
    else max_eta_h = tauh_eta;
    if (VBFjet1_eta > VBFjet2_eta) max_eta_vbf = VBFjet1_eta; 
    else max_eta_vbf = VBFjet2_eta;
    return max_eta_vbf - max_eta_h;
}


inline float FeatComp::calc_mt_tot(const LorentzVector& l_1, const LorentzVector& l_2, const LorentzVector& met) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */

    return std::sqrt(std::pow(FeatComp::calc_mt(l_1, met), 2) + std::pow(FeatComp::calc_mt(l_2, met), 2) + std::pow(FeatComp::calc_mt(l_1, l_2), 2));
}


float FeatComp::calc_pzeta(const LorentzVector& l_1, const LorentzVector& l_2, const LorentzVector& met) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */

    const LorentzVector ll_p4 = l_1 + l_2;
    const TVector2 ll_p2(ll_p4.Px(), ll_p4.Py());
    const TVector2 met_p2(met.Px(), met.Py());
    const TVector2 ll_s = ll_p2 + met_p2;
    const TVector2 l1_u(std::cos(l_1.Phi()), std::sin(l_1.Phi()));
    const TVector2 l2_u(std::cos(l_2.Phi()), std::sin(l_2.Phi()));
    const TVector2 ll_u = l1_u + l2_u;
    const float ll_u_met = ll_s * ll_u;
    const float ll_mod = ll_u.Mod();
    return ll_u_met / ll_mod;
}


float FeatComp::calc_pzeta_visible(const LorentzVector& l_1, const LorentzVector& l_2) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */

    const LorentzVector ll_p4 = l_1 + l_2;
    const TVector2 ll_p2(ll_p4.Px(), ll_p4.Py());
    const TVector2 l1_u(std::cos(l_1.Phi()), std::sin(l_1.Phi()));
    const TVector2 l2_u(std::cos(l_2.Phi()), std::sin(l_2.Phi()));
    const TVector2 ll_u = l1_u + l2_u;
    const float ll_p2u = ll_p2 * ll_u;
    const float ll_mod = ll_u.Mod();
    return ll_p2u / ll_mod;
}


std::pair<float, float> FeatComp::calc_top_masses(const LorentzVector& l_1, const LorentzVector& l_2, const LorentzVector& b_1, const LorentzVector& b_2,
                                                  const LorentzVector& met) {
    /* Modified from https://github.com/hh-italian-group/AnalysisTools/blob/1be0da0748d69827ed7ebda6d9b8198b87f170fd/Core/include/AnalysisMath.h */

    std::vector<std::pair<float, float>> vector_mass_top = {
        {(l_1 + b_1 + met).mass(), (l_2 + b_2).mass()},
        {(l_1 + b_2 + met).mass(), (l_2 + b_1).mass()},
        {(l_1 + b_1).mass(),       (l_2 + b_2 + met).mass()},
        {(l_1 + b_2).mass(),       (l_2 + b_1 + met).mass()}
    };
    std::vector<std::pair<unsigned int, float>> distance;
    for (unsigned int i = 0; i < vector_mass_top.size(); ++i) distance.emplace_back(i, pow(vector_mass_top[i].first - 172.5,2) + pow(vector_mass_top[i].second - 172.5,2));
    std::sort(distance.begin(), distance.end(), [](const std::pair<unsigned int, float>& el1,const std::pair<unsigned int, float>& el2) {
        return el1.second < el2.second;
    });
    return vector_mass_top.at(distance.front().first);
}