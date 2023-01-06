@quickactivate "RasterUniqueCharacterizationCode"

@warn "TRIVIAL TEST PARAMS"
n_trials=90;
n_tests=3;
n_signals=13;
MEA_spikes=167;
lag_extents=(10,10);
signal_jitters=(12,12);

include(scriptsdir("detection","trial_detect.jl"))