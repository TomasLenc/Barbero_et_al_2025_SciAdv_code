function [s_all, ratios, iois] = make_slice_signals(fs, event, n_cycles, cycle_dur, varargin)

parser = inputParser; 

addParameter(parser, 'vol', []); 
addParameter(parser, 'ratios', []); 
addParameter(parser, 'proto1', []); 
addParameter(parser, 'proto2', []); 
addParameter(parser, 'n_cond', []); 
addParameter(parser, 'n_cond_beyond_proto1', 0); 
addParameter(parser, 'n_cond_beyond_proto2', 0); 

parse(parser, varargin{:}); 

vol = parser.Results.vol; 
ratios = parser.Results.ratios; 
proto1 = parser.Results.proto1; 
proto2 = parser.Results.proto2; 
n_cond = parser.Results.n_cond; 
n_cond_beyond_proto1 = parser.Results.n_cond_beyond_proto1; 
n_cond_beyond_proto2 = parser.Results.n_cond_beyond_proto2; 

%%

if isempty(vol)
    vol = ones(size(proto1)); 
end
if isempty(ratios)
    ratios = interpolate_prototypes(proto1, proto2, n_cond, cycle_dur, ...
                                    'n_cond_beyond_proto1', n_cond_beyond_proto1, ...
                                    'n_cond_beyond_proto2', n_cond_beyond_proto2); 
end

n_cond_total = size(ratios, 1);

s_all = cell(1, n_cond_total); 

iois = ratios * cycle_dur; 

N = round(n_cycles * cycle_dur * fs); 

% Loop over conditions (interval ratio)
for i_cond = 1:n_cond_total

    s = zeros(1, N); 
    
    for i_cycle = 1:n_cycles

        t_onsets = iois_to_onset_times(iois(i_cond, :), ...
                                       cycle_dur, ...
                                       i_cycle); 

        idx_onsets = round(t_onsets * fs); 

        s = insert_sound_events(s, idx_onsets, event, 'vol', vol);

    end

    s_all{i_cond} = s(1:N); 
    
end
    

