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
#include <Math/PxPyPzM4D.h>

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
										 const float& c3,
										 const bool& cut_pass);
	std::vector<float> process_as_vec(const LorentzVector& b_1,
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
									  const float& c3,
									  const bool& cut_pass);
	void process_to_vec(std::vector<std::unique_ptr<float>>& feats,
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
						const float& c3,
						const bool& cut_pass);
	std::vector<std::string> get_feats();
};

#endif /* EVT_PROC_HH_ */
