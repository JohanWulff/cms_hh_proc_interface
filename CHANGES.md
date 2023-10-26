# V4.0

- `EvtProc::process`, `EvtProc::process_as_vec`, and `EvtProc::process_to_vec` now take a `cut_pass` Boolean argument indicating is the event passes the usual ellipse/mass-window cut

# V3.0

- `FeatComp::process` now takes the following extra arguments:
	- const float& b_1_cvsl,
	- const float& b_2_cvsl,
	- const float& vbf_1_cvsl,
	- const float& vbf_2_cvsl,
	- const float& b_1_cvsb,
	- const float& b_2_cvsb,
	- const float& vbf_1_cvsb,
	- const float& vbf_2_cvsb
- Fixed working points for all taggers, per year
- Charm tagging features now categorical rather than float, float values are available as `*_raw` for feature selection
- `HHKin_mass_raw` now set to zero rather than NaN when invalid

# V2.0

- `FeatComp::process`, `EvtProc::process`, `EvtProc::process_as_vec`, and `EvtProc::process_to_vec` now only take a single set of CSV values, rather than CSV and deep*. Pass the CSV values you wish to use. 
- `FeatComp::FeatComp` and `EvtProc::EvtProc` `_use_deep_csv` arguments renamed to `_use_deep_bjet_wps`
- `EvtProc::process`, `EvtProc::process_as_vec`, and `EvtProc::process_to_vec` now take the following extra arguments:
    - const float& b_1_hhbtag,
	- const float& b_2_hhbtag,
	- const float& vbf_1_hhbtag,
	- const float& vbf_2_hhbtag,
    - const float& b_1_cvsl,
	- const float& b_2_cvsl,
	- const float& vbf_1_cvsl,
	- const float& vbf_2_cvsl,
	- const float& b_1_cvsb,
	- const float& b_2_cvsb,
	- const float& vbf_1_cvsb,
	- const float& vbf_2_cvsb,
    - const float& cv,
	- const float& c2v,
	- const float& c3
- `EvtProc::process`, `EvtProc::process_as_vec`, and `EvtProc::process_to_vec` now no longer take the following arguments:
	- const float& mt_tot,
    - const float& p_zetavisible,
    - const float& p_zeta,
    - const float& top_1_mass,
    - const float& top_2_mass,
    - const float& l_1_mt,
    - const float& l_2_mt,
