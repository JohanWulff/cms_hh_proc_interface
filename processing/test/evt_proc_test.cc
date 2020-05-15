#include "cms_hh_proc_interface/processing/interface/feat_comp.hh"
#include "cms_hh_proc_interface/processing/interface/evt_proc.hh"
#include <iostream>
#include <string>
#include <map>
#include <random>
#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PxPyPzM4D.h>

int main(int argc, char *argv[]) {
	using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> mom(-100.,100.);
    std::uniform_real_distribution<> energy(20.,500.);
    std::uniform_real_distribution<> csv(-1.,1.);

    std::cout << "Generating random event... ";
    LorentzVector b_1(mom(rng),mom(rng),mom(rng),energy(rng));
    LorentzVector b_2(mom(rng),mom(rng),mom(rng),energy(rng));
    LorentzVector l_1(mom(rng),mom(rng),mom(rng),energy(rng));
    LorentzVector l_2(mom(rng),mom(rng),mom(rng),energy(rng));
    LorentzVector met(mom(rng),mom(rng),0,       energy(rng));
    LorentzVector sv(mom(rng),mom(rng),mom(rng),energy(rng));
    LorentzVector vbf_1(mom(rng),mom(rng),mom(rng),energy(rng));
    LorentzVector vbf_2(mom(rng),mom(rng),mom(rng),energy(rng));
    float hh_kinfit_mass = energy(rng);
    float hh_kinfit_chi2 = energy(rng);
    float mt2            = energy(rng);
    float mt_tot         = energy(rng);
    float pzetavisible   = energy(rng);
    float pzeta          = energy(rng);
    float top_1_mass     = energy(rng);
    float top_2_mass     = energy(rng);
    float l_1_mt         = energy(rng);
    float l_2_mt         = energy(rng);
    bool is_boosted      = csv(rng) > 1.;
    float csv_1(csv(rng)), csv_2(csv(rng));
    float b_1_hhbtag(csv(rng)), b_2_hhbtag(csv(rng)), vbf_1_hhbtag(csv(rng)), vbf_2_hhbtag(csv(rng));
    Channel channel = tauTau;
    Year year = y16;
    float res_mass = 450;
    Spin spin = nonres;
    float klambda = 1;
    std::cout << "Generated\n";

    std::cout << "Instantiating return_all EvtProc... ";
    EvtProc evt_proc;
    std::cout << "Instantiated\n ";
    std::cout << "Processing event... ";
    std::map<std::string, float> feats = evt_proc.process(b_1, b_2, l_1, l_2, met, sv, vbf_1, vbf_2, hh_kinfit_mass, hh_kinfit_chi2, mt2, mt_tot, pzetavisible,
                                                          pzeta, top_1_mass, top_2_mass, l_1_mt, l_2_mt, is_boosted, csv_1, csv_2,
                                                          channel, year, res_mass, spin, klambda, 1, true, false,
                                                          b_1_hhbtag, b_2_hhbtag, vbf_1_hhbtag, vbf_2_hhbtag);
    std::cout << "Processed\n";
    for (auto const& f : feats) std::cout << f.first << " : " << f.second << "\n";
    std::cout << feats.size() << " features returned\n\n";

    std::cout << "Instantiating requested EvtProc... ";
    std::vector<std::string> requested = {"costheta_met_htt", "phi", "hh_kinfit_chi2", "boosted", "costheta_l1_httmet"};
    EvtProc evt_proc2(false, requested);

    std::cout << "\nProcessing event as vector... ";
    std::vector<float> vec = evt_proc2.process_as_vec(b_1, b_2, l_1, l_2, met, sv, vbf_1, vbf_2, hh_kinfit_mass, hh_kinfit_chi2, mt2, mt_tot, pzetavisible,
                                                      pzeta, top_1_mass, top_2_mass, l_1_mt, l_2_mt, is_boosted, csv_1, csv_2, channel,
                                                      year, res_mass, spin, klambda, 2, false, true, b_1_hhbtag, b_2_hhbtag, vbf_1_hhbtag, vbf_2_hhbtag);
    std::cout << "Processed\n";
    std::cout << "Recieved vector of " << vec.size() << " elements\n";
    for (auto const& f : vec) std::cout << f << " ";
    std::cout << "\n";

    std::cout << "Checking feat names... ";
    std::vector<std::string> names = evt_proc2.get_feats();
    std::cout << "Processed\n";
    std::cout << "Recieved:\n";
    for (auto const& f : names) std::cout << f << " ";
    std::cout << "\n";

    std::cout << "Processing event to vector of pointers to zeros... ";
    unsigned int _n_feats = names.size();
    std::vector<std::unique_ptr<float>> feat_vals;
    feat_vals.reserve(_n_feats);
    for (unsigned int i = 0; i < _n_feats; i++) feat_vals.emplace_back(new float(0));
    evt_proc2.process_to_vec(feat_vals, b_1, b_2, l_1, l_2, met, sv, vbf_1, vbf_2, hh_kinfit_mass, hh_kinfit_chi2, mt2, mt_tot, pzetavisible, pzeta,
                             top_1_mass, top_2_mass, l_1_mt, l_2_mt, is_boosted, csv_1, csv_2, channel, year, res_mass, spin, klambda,
                             3, false, false, b_1_hhbtag, b_2_hhbtag, vbf_1_hhbtag, vbf_2_hhbtag);
    std::cout << "Processed, vector values now: \n";
    for (unsigned int i = 0; i < _n_feats; i++) std::cout << *(feat_vals[i]) << " ";
    std::cout << "\n";
    
    return 0;
}