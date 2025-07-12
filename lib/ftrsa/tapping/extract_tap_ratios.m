function [tap_ratios, tap_ratios_mean, n_valid_tap_ratios, n_possible_tap_ratios, tap_ratios_all] = ...
                        extract_tap_ratios(tap_onsets, stim_iois, n_cycles, varargin)
% Extract tap times from continuous time series (e.g. recorded using a
% microphone). The extraction is based on thresholding the signal's
% amplitude. 
% 
% Parameters
% ----------
% x : shape=[1, time]
%     Input time series. 
% fs : int
%     Sampling rate.
% stim_iois : array_like, shape=[1, N_iois_per_cycle]
%     Assuming the tapping was a response to a stimulus sequence made of a
%     repeating rhythmic pattern, this vector describes the
%     inter-onset intervals between successive sound comprising one
%     repeating pattern in the stimulus sequence. First element is the time
%     between the sound starting the repeating rhythmic cycle and the
%     following sound, the second element is the time to the next sound
%     etc. etc. The function can deal with any N of elements within the
%     rhythmic cycle. 
% tap_tone_min_asy : float, optional, default=$
%     A tap onset is assigned to a particular stimulus sound event if it is
%     less than `tap_tone_min_asy` away from that event. 
% n_cycles_skip : int, optional, default=$
%     Number of cycles at the begining of the sequence that will be
%     ignorede. This can be used if we want to discard "unstable" tapping
%     at the begining of the trial. 
% plot_diagnostic : bool, optional, default=false
%     If true, some diagnostic figures will be generated when the function
%     is called, so that you can visually check whether the onset
%     extraction was reasonably good given the threshold parameters etc. 
% verbose : bool, optional, default=true
%     Do you want many warnings? :) 
% 
% Returns 
% -------
% tap_ratios : array_like, shape=[N_valid_cycles, N_iois_per_cycle]
%     Ratios observed in the tapping response. Only ratios for valid cycles
%     are returned. The response to a cycle, i.e. one repetition of the
%     rhytjhmic pattern constituting the stimulus sequence, is considered
%     "valid" if there is a tap onset assigned to each sound event in that
%     cycle AND the starting event marking the onset of the following
%     cycle. 
% tap_ratios_mean : array_like, shape=[1, N_iois_per_cycle]
%     Average of tap_ratios. 
% n_valid_tap_ratios : int
%     Number of cycles where the tapping is considered "valid". 
% n_possible_tap_ratios : int
%     Maximum numnber of cycles where the tapping could considered "valid". 
% tap_ratios_all : array_like, shape=[N_cycles, N_iois_per_cycle]
%     Ratios observed in the tapping response for all cycles. If the tap
%     ratio in the particular cycle was not valid, NaN will be inserted in
%     the corresponding row.
% 

parser = inputParser; 

addParameter(parser, 'tap_tone_min_asy', 0.050); 
addParameter(parser, 'mean_asy', 0); 
addParameter(parser, 'n_cycles_skip', 0); 
addParameter(parser, 'plot_diagnostic', false); 
addParameter(parser, 'verbose', true); 

parse(parser, varargin{:}); 

tap_tone_min_asy = parser.Results.tap_tone_min_asy; 
mean_asy = parser.Results.mean_asy; 
n_cycles_skip = parser.Results.n_cycles_skip; 
plot_diagnostic = parser.Results.plot_diagnostic; 
verbose = parser.Results.verbose; 

%%

cycle_dur = sum(stim_iois); 

sound_onsets0 = [cycle_dur * n_cycles_skip : cycle_dur : cycle_dur * (n_cycles - 1)]; 

sound_onsets = nan(length(sound_onsets0), length(stim_iois) + 1); 

sound_onsets(:, 1) = sound_onsets0; 

for i_onset=2:size(sound_onsets, 2)
    sound_onsets(:, i_onset) = sound_onsets0 + sum(stim_iois(1 : (i_onset-1)));     
end


tap_ratios_all = nan(size(sound_onsets, 1), 2); 

% go over pattern cycles 
for i_cycle=1:size(sound_onsets, 1)

    % find the closest tap to the start of the cycle and each sound onset
    asy = nan(1, size(sound_onsets, 2)); 
    asy_adj = nan(1, size(sound_onsets, 2)); 
    idx = nan(1, size(sound_onsets, 2)); 
    
    % !! we can in principle assign the same tap to two sound onsets....
    % this can be fixed in the future: greedy style 
    for i_onset=1:size(sound_onsets, 2)
        
        % account for mean asynchrony when picking the closest tap
        [~, idx(i_onset)] = ...
                    min(abs((tap_onsets-mean_asy) - sound_onsets(i_cycle, i_onset))); 
                
        asy(i_onset) = tap_onsets(idx(i_onset)) - sound_onsets(i_cycle, i_onset); 
        asy_adj(i_onset) = asy(i_onset) - mean_asy; 
        
    end
        
    % if there are valid taps for each of the tones, accept this
    % cycle and compute the ratio 
    if all( abs(asy_adj) < tap_tone_min_asy )

        if any((idx + 1) > length(tap_onsets)) 
            if verbose
                warning('cycle %d last tap selected as a start of interval...skipping', i_cycle); 
            
            end
            continue
        end
        
        dur_cycle = tap_onsets(idx(end)) - tap_onsets(idx(1)); 
        
        dur = nan(1, length(stim_iois)); 
        
        for i_iti=1:length(stim_iois)
            dur(i_iti) = tap_onsets(idx(i_iti)+1) - tap_onsets(idx(i_iti));         
        end
        
        tap_ratios_all(i_cycle, :) = dur ./ dur_cycle;
    end

end

% get only valid tap ratios
tap_ratios = tap_ratios_all(~isnan(tap_ratios_all(:,1)), :); 

tap_ratios_mean = mean(tap_ratios, 1); 
n_valid_tap_ratios = size(tap_ratios, 1);     
n_possible_tap_ratios = size(sound_onsets, 1); 

%%

% sanity check
if plot_diagnostic

    figure('color', 'w', 'pos', [119 682 1662 126]); 
    plot([sound_onsets(:, 1), sound_onsets(:, 1)], [-1, 1], ':k')
    hold on 
    for i=2:size(sound_onsets, 2)
        plot([sound_onsets(:, i), sound_onsets(:, i)], [-1, 1], ':b')
    end
    plot([tap_onsets], [0], 'ro')
    box off
    ylim([-1, 1])
end

