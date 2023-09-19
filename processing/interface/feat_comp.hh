#ifndef FEAT_COMP_HH_
#define FEAT_COMP_HH_

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
#include <Math/PxPyPzM4D.h>
#include <Math/PtEtaPhiE4D.h>
#include <TMath.h>
#include <TVector2.h>

enum Channel{tauTau=0, muTau=1, eTau=2};
enum Year{y16=0, y17=1, y18=2, y16APV=3};

class FeatComp {
	/* Class for computing requested features for final-state LorentzVectors */

private:
	// Names
	using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

	// Variables
    bool _all, _use_deep_bjet_wps;
    std::vector<std::string> _requested;
	std::map<int,std::vector<float>> _bjet_wps{{0, {0,0.5803,0.8838,0.9693}},   // Not provided for 2016, copied 2017
											   {1, {0,0.5803,0.8838,0.9693}},
											   {2, {0,0.5803,0.8838,0.9693}}};  // Not provided for 2018, copied 2017
	std::map<int,std::vector<float>> _deep_bjet_wps{{0, {0,0.0614,0.3093,0.7221}},
													{1, {0,0.0532,0.4506,0.7738}},
													{2, {0,0.0490,0.2783,0.7100}}};
	std::map<int,std::vector<float>> _cvsl_wps{{0, {0,0.03,0.085,0.48}},
											   {1, {0,0.03,0.085,0.520}},  // Not provided for 2017, copied 2016
											   {2, {0,0.038,0.099,0.282}}};	
	std::map<int,std::vector<float>> _cvsb_wps{{0, {0.4,0.29,0.05,0}},
											   {1, {0.4,0.34,0.05,0}},
											   {2, {0.246,0.325,0.267,0}}};

	// Methods
	void _add_btag_flags(Year year, const float& b_1_csv, const float& b_2_csv, std::map<std::string, float>& feats);
	int _get_cvsl_flag(Year year, const float& score);
	int _get_cvsb_flag(Year year, const float& score);
	inline bool _feat_check(std::string);

public:
    // Methods
	FeatComp(bool return_all=true, std::vector<std::string> requested={}, bool use_deep_bjet_wps=true);
	~FeatComp();
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
                                               const float& DeepMET_ResolutionTune_py);
	inline float delta_eta(const LorentzVector&, const LorentzVector&);
	inline float delta_phi(const LorentzVector&, const LorentzVector&);
	inline float delta_r(const LorentzVector&, const LorentzVector&);
	inline float delta_r_boosted(const LorentzVector&, const LorentzVector&, const LorentzVector&);
	inline float calc_mt(const LorentzVector&, const LorentzVector&);
	inline float calc_phi(const LorentzVector&, const LorentzVector&, const LorentzVector&, const LorentzVector&, const LorentzVector&);
	inline float calc_phi_1(const LorentzVector&, const LorentzVector&, const LorentzVector&, const LorentzVector&);
	inline float calc_cos_delta_star(const LorentzVector&, const LorentzVector&);
	inline float calc_cos_delta(const LorentzVector&, const LorentzVector&);
  	inline float calc_centrality( const LorentzVector&, const LorentzVector&, const LorentzVector&);
  	inline float calc_hh_centrality( const LorentzVector&, const LorentzVector&, const LorentzVector&, const LorentzVector&);
  	inline float calcDeltaEtaMinus( const LorentzVector&, const LorentzVector&, const LorentzVector&, const LorentzVector&);
  	inline float calcDeltaEtaPlus( const LorentzVector&, const LorentzVector&, const LorentzVector&, const LorentzVector&);
  	inline float calc_mt_tot(const LorentzVector& l_1, const LorentzVector& l_2, const LorentzVector& met);
	inline float calc_dmet_px(const float& px,const float& py,const float& phi);
	inline float calc_dmet_py(const float& px,const float& py,const float& phi);
	float mpi_to_pi(const float& phi);
	float calc_pzeta(const LorentzVector& l_1, const LorentzVector& l_2, const LorentzVector& met);
	float calc_pzeta_visible(const LorentzVector& l_1, const LorentzVector& l_2);
	std::pair<float, float> calc_top_masses(const LorentzVector& l_1, const LorentzVector& l_2, const LorentzVector& b_1, const LorentzVector& b_2,
                                          const LorentzVector& met);
};

#endif /* FEAT_COMP_HH_ */
