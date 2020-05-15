# May

- `FeatComp::process`, `EvtProc::process`, `EvtProc::process_as_vec`, and `EvtProc::process_to_vec` now only take a single set of CSV values, rather than CSV and deep*. Pass the CSV values you wish to use. 
- `FeatComp::FeatComp` and `EvtProc::EvtProc` `_use_deep_csv` arguments renamed to `_use_deep_bjet_wps`
