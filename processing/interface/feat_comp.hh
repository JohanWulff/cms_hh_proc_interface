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
#include <TVector2.h>

enum Channel{tauTau=0, muTau=1, eTau=2};
enum Year{y16=0, y17=1, y18=2};

class FeatComp {
	/* Class for computing requested features for final-state LorentzVectors */

private:
	// Names
	using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

	// Variables
    bool _all, _use_deep_bjet_wps;
    std::vector<std::string> _requested;
	std::vector<float> _bjet_wps{0, 0.5426, 0.8484, 0.9535};  // TODO: Check these
	std::vector<float> _deep_bjet_wps{0, 0.0494, 0.2770, 0.7264};  //TODO: Update these

	// Methods
	void _add_jet_flags(const float& b_1_csv, const float& b_2_csv, std::map<std::string, float>& feats);
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
										 const float& hh_kinfit_m,
                                         const bool& is_boosted,
                                         const float& b_1_csv,
                                         const float& b_2_csv,
                                         Channel channel,
                                         Year year,
                                         const int& n_vbf,
                                         const bool& svfit_conv,
                                         const bool& hh_kinfit_conv);
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
  float calc_pzeta(const LorentzVector& l_1, const LorentzVector& l_2, const LorentzVector& met);
  float calc_pzeta_visible(const LorentzVector& l_1, const LorentzVector& l_2);
  std::pair<float, float> calc_top_masses(const LorentzVector& l_1, const LorentzVector& l_2, const LorentzVector& b_1, const LorentzVector& b_2,
                                          const LorentzVector& met);
};

#endif /* FEAT_COMP_HH_ */
