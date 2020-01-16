#include "cms_hh_proc_interface/processing/interface/feat_comp.hh"
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

    std::cout << "Instantiating FeatComp... ";
    FeatComp feat_comp;
    std::cout << "Instantiated\n ";

    std::cout << "Generating random event... ";
    LorentzVector b_1(mom(rng),mom(rng),mom(rng),energy(rng));
    LorentzVector b_2(mom(rng),mom(rng),mom(rng),energy(rng));
    LorentzVector l_1(mom(rng),mom(rng),mom(rng),energy(rng));
    LorentzVector l_2(mom(rng),mom(rng),mom(rng),energy(rng));
    LorentzVector met(mom(rng),mom(rng),0,       energy(rng));
    LorentzVector sv(mom(rng),mom(rng),mom(rng),energy(rng));
    float hh_kinfit_mass = energy(rng);
    bool is_boosted = csv(rng) > 1.;
    float csv_1(csv(rng)), csv_2(csv(rng)), deepcsv_1(csv(rng)), deepcsv_2(csv(rng));
    Channel channel = tauTau;
    std::cout << "Generated\n";

    std::cout << "Processing event... ";
    std::map<std::string, float> feats = feat_comp.process(b_1, b_2, l_1, l_2, met, sv, hh_kinfit_mass, is_boosted,
                                                           csv_1, csv_2, deepcsv_1, deepcsv_2, channel);
    std::cout << "Processed\n";

    for (auto const& f : feats) std::cout << f.first << ":" << f.second << "\n";
}