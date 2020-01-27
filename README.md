# cms_hh_proc_interface
Interface providing feature computation from ROOT LorentzVectors for CMS RunII hh->bbtautau DNN

# Installation

1. cmsrel CMSSW_10_2_15
1. cd CMSSW_10_2_15/src
1. cmsenv
1. git clone git@github.com:GilesStrong/cms_hh_proc_interface.git
1. scram b -j 12

# Testing

Building the tool in CMSSW will generate two executables in `CMSSW_10_2_16/test/slc7_amd64_gcc700/` called: `testfeatcomp` which will generate the example and test the feature computation object; and `testevtproc` which will check the event processing object.

# Usage

The interface consists of two objects:
- `FeatComp` computes features and categorical values from `ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>` objects and other event information and returns it as a map of feature name to float value (`FeatComp::process`). If `all` is set to `true` then all features will be computed, otherwise a set of feature names can be passed to `requested`, in which case only the features listed will be computed (**NB** the order of the returned map is not guaranteed to match the order of requested features!). The choice of whether to use CSV or Deep CSV can can be controlled using `use_deep_csv`.
- `EvtProc` is the main interface for processing events from `ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>` objects and other event information by wrapping `FeatComp` and adding additional information (`EvtFeat::process`). If `all` is set to `true` then all features will be computed, otherwise a set of feature names can be passed to `requested`, in which case only the features listed will be computed and **returned in the requested order**. The choice of whether to use CSV or Deep CSV can can be controlled using `use_deep_csv`. If a vector of floats is required, i.e. for passing straight to a DNN, then `EvtFeat::process_as_vec` can be called instead.

# Notes

- The enum objects `Channel` and `Year` are defined in `feat_comp.hh`.
