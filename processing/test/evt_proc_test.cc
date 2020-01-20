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
    float csv_1(csv(rng)), csv_2(csv(rng)), deepcsv_1(csv(rng)), deepcsv_2(csv(rng));
    Channel channel = tauTau;
    Year year = y16;
    float res_mass = 450;
    std::cout << "Generated\n";

    std::cout << "Instantiating return_all EvtProc... ";
    EvtProc evt_proc;
    std::cout << "Instantiated\n ";

    std::cout << "Processing event... ";
    std::map<std::string, float> feats = evt_proc.process(b_1, b_2, l_1, l_2, met, sv, hh_kinfit_mass, hh_kinfit_chi2, mt2, mt_tot, pzetavisible, pzeta,
                                                          top_1_mass, top_2_mass, l_1_mt, l_2_mt,
                                                          is_boosted, csv_1, csv_2, deepcsv_1, deepcsv_2, channel, year, res_mass);
    std::cout << "Processed\n";
    for (auto const& f : feats) std::cout << f.first << " : " << f.second << "\n";
    std::cout << feats.size() << " features returned\n\n";

    std::cout << "Instantiating requested EvtProc... ";
    std::set<std::string> requested = {"costheta_met_htt", "phi", "hh_kinfit_chi2", "is_boosted", "costheta_l1_httmet"};
    EvtProc evt_proc2(false, requested);
    std::cout << "Instantiated\n ";

    std::cout << "Processing event... ";
    feats = evt_proc2.process(b_1, b_2, l_1, l_2, met, sv, hh_kinfit_mass, hh_kinfit_chi2, mt2, mt_tot, pzetavisible, pzeta,
                              top_1_mass, top_2_mass, l_1_mt, l_2_mt,
                              is_boosted, csv_1, csv_2, deepcsv_1, deepcsv_2, channel, year, res_mass);
    std::cout << "Processed\n";

    std::cout << "Expected: ";
    for (auto const& f : requested) std::cout << f << " ";
    std::cout << "\nReceived: ";
    for (auto const& f : feats)     std::cout << f.first << " ";

    return 0;
}