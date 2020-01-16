#ifndef FEAT_COMP_HH_
#define FEAT_COMP_HH_

// C++
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>

// ROOT
#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PxPyPzM4D.h>

using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

enum Channel{tauTau=0, muTau=1, eTau=2};

class FeatComp {
	/* Class for computing requested features for final-state LorentzVectors*/

private:
	// Variables
    bool _all, _verbose, _use_deep_csv;
    std::set<std::string> _requested;
	std::vector<float> _csv_wps{0, 0.5426, 0.8484, 0.9535};  // TODO: Check these
	std::vector<float> _deep_csv_wps{0, 0.5426, 0.8484, 0.9535};  //TODO: Update these

	// Methods
	void _add_jet_flags(const float&, const float&, const float&, const float&, std::map<std::string, float>&);
	inline bool _feat_check(std::string);

public:
    // Methods
	FeatComp(bool return_all=true, std::set<std::string> requested={}, bool use_deep_csv=true, bool verbose=false);
	~FeatComp();
	std::map<std::string, float> process(const LorentzVector&,  // b_1
										 const LorentzVector&,  // b_2
										 const LorentzVector&,  // l_1
										 const LorentzVector&,  // l_2
										 const LorentzVector&,  // MET
										 const LorentzVector&,  // SVFit
										 const float&,          // HH KinFit mass
										 const bool&, 			// Is boosted
										 const float&,			// b_1 CSV
										 const float&,			// b_2 CSV
										 const float&,			// b_1 Deep CSV
										 const float&,			// b_2 Deep CSV
										 Channel);	            // Channel		
	inline float delta_eta(const LorentzVector&, const LorentzVector&);
	inline float delta_phi(const LorentzVector&, const LorentzVector&);
	inline float delta_r(const LorentzVector&, const LorentzVector&);
	inline float delta_r_boosted(const LorentzVector&, const LorentzVector&, const LorentzVector&);
	inline float calc_mt(const LorentzVector&, const LorentzVector&);
	inline float calc_phi(const LorentzVector&, const LorentzVector&, const LorentzVector&, const LorentzVector&, const LorentzVector&);
	inline float calc_phi_1(const LorentzVector&, const LorentzVector&, const LorentzVector&, const LorentzVector&);
	inline float calc_cos_delta_star(const LorentzVector&, const LorentzVector&);
	inline float calc_cos_delta(const LorentzVector&, const LorentzVector&);
};

#endif /* FEAT_COMP_HH_ */
