function [asy] = get_asynchronies(tap_onsets, stim_iois, n_cycles, varargin)
% For each tap, find the closest sound onset time and compute signed
% asynchrony from it. 

parser = inputParser; 

addParameter(parser, 'tap_tone_min_asy', 0.050); 
addParameter(parser, 'n_cycles_skip', 0); 
addParameter(parser, 'plot_diagnostic', false); 
addParameter(parser, 'verbose', true); 

parse(parser, varargin{:}); 

n_cycles_skip = parser.Results.n_cycles_skip; 

%%

cycle_dur = sum(stim_iois); 

sound_onsets = repmat([0, stim_iois(1)], n_cycles, 1) + cycle_dur * [n_cycles_skip : n_cycles-1]'; 

sound_onsets = reshape(sound_onsets', [], 1); 

% find the closest tap to the start of the cycle and each sound onset
asy = nan(1,length(tap_onsets)); 

for i_tap=1:length(tap_onsets)

    [~, idx_closest_sound] = min(abs(tap_onsets(i_tap) - sound_onsets)); 
    
    asy(i_tap) = tap_onsets(i_tap) - sound_onsets(idx_closest_sound); 
            
end





