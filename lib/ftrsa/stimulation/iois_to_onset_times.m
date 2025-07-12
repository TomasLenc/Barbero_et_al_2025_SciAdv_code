function t_onsets = iois_to_onset_times(iois, pattern_period, i_cycle)

t_from_cycle_onset = cumsum(iois); 

t_onsets = (i_cycle-1) * pattern_period + [0, t_from_cycle_onset(1:end-1)]; 
