#ifndef EVT_PROC_HH_
#define EVT_PROC_HH_

// C++
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <algorithm>

// ROOT
#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PxPyPzE4D.h>
#include <Math/PtEtaPhiE4D.h>

// Local
#include "feat_comp.hh"

enum Spin{radion=0, graviton=1, nonres=2};

class EvtProc {
	/* Class for extracting and computing required information and returning it expected order */

private:
    // Names
    using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

	// Variables
    bool _all;
	FeatComp* _feat_comp;
    std::vector<std::string> _requested;

	// Methods
	inline bool _feat_check(std::string);

public:
    // Methods
	EvtProc(bool return_all=true, std::vector<std::string> requested={}, bool use_deep_bjet_wps=true);
	~EvtProc();
	std::map<std::string, float> process(const LorentzVector& b_1,
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
                                              const bool&  isBoosted,
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
                                              const float& bjet2_cID);
	std::vector<float> process_as_vec(const LorentzVector& b_1,
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
                                              const bool&  isBoosted,
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
                                              const float& bjet2_cID);
	void process_to_vec(std::vector<std::unique_ptr<float>>& feats,
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
                                              const bool&  isBoosted,
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
                                              const float& bjet2_cID);
	std::vector<std::string> get_feats();
};

#endif /* EVT_PROC_HH_ */
