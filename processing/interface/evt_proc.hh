#ifndef EVT_PROC_HH_
#define EVT_PROC_HH_

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

// Local
#include "feat_comp.hh"

class EvtProc {
	/* Class for extracting and computing required information and returning it expected order */

private:
    // Names
    using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

	// Variables
    bool _all;
	FeatComp* _feat_comp;
    std::set<std::string> _requested;

	// Methods
	inline bool _feat_check(std::string);
	std::map<std::string, float> _sort_feats(std::map<std::string, float>);

public:
    // Methods
	EvtProc(bool return_all=true, std::set<std::string> requested={}, bool use_deep_csv=true);
	~EvtProc();
	std::map<std::string, float> process(const LorentzVector&,  // b_1
										 const LorentzVector&,  // b_2
										 const LorentzVector&,  // l_1
										 const LorentzVector&,  // l_2
										 const LorentzVector&,  // MET
										 const LorentzVector&,  // SVFit
										 const float&,          // HH KinFit mass
										 const float&,          // HH KinFit chi2
										 const float&,          // MT2
										 const float&,          // MT total
										 const bool&, 			// Is boosted
										 const float&,			// b_1 CSV
										 const float&,			// b_2 CSV
										 const float&,			// b_1 Deep CSV
										 const float&,			// b_2 Deep CSV
										 Channel,	            // Channel
										 const float&);		    // Resonant mass (set to zero is non-resonant)
};

#endif /* EVT_PROC_HH_ */
