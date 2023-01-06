@quickactivate "RasterUniqueCharacterizationCode"

n_trials=15;
n_tests=15;
n_signals=13;
MEA_spikes=167;
lag_extents=(10,10);
signal_jitters=(12,12);

include(scriptsdir("detection","trial_detect_RATE.jl"))
zscore!