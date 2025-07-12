function [ratios, iois] = interpolate_prototypes(proto1, proto2, n_cond, pattern_period, varargin)

parser = inputParser; 

addParameter(parser, 'n_cond_beyond_proto1', 0); 
addParameter(parser, 'n_cond_beyond_proto2', 0); 
addParameter(parser, 'interp_type', 'line'); % line or arch 
addParameter(parser, 'direction', 'anticlockwise'); % only applies to arch 

parse(parser, varargin{:}); 

n_cond_beyond_proto1 = parser.Results.n_cond_beyond_proto1; 
n_cond_beyond_proto2 = parser.Results.n_cond_beyond_proto2; 
interp_type = parser.Results.interp_type; 
direction = parser.Results.direction; 

%%

if strcmp(interp_type, 'line')
    n_intervals = length(proto1); 

    proto1 = proto1 ./ sum(proto1); 
    proto2 = proto2 ./ sum(proto2); 

    diff_vec = proto2 - proto1; 

    diff_vec_to_add = diff_vec / (n_cond-1); 

    n_diff_to_add = [-n_cond_beyond_proto1 : (n_cond - 1 + n_cond_beyond_proto2)]; 

    ratios = nan(length(n_diff_to_add) , n_intervals); 
    iois = nan(length(n_diff_to_add) , n_intervals); 

    for iC=1:length(n_diff_to_add) 

        % compute ratios between the 3 time intervals in the pattern
        ratios(iC, :) = proto1 + n_diff_to_add(iC) * diff_vec_to_add; 

        % convert ratios to seconds
        iois(iC, :) = pattern_period * ratios(iC, :); 

        % sanity check 
        assert(sum(iois(iC, :)) - pattern_period < 1e12); 
    end


elseif strcmp(interp_type, 'arch')

    % !!! this ONLY works for 3ioi rhythms!!!
    if numel(proto1) ~=3 ||  numel(proto2) ~=3
        error('Requested arch interpolation but this requires 3-interval rhyrhms!!!'); 
    end
    
    % find angle pointing from 1,1,1 to proto1 
    q = tri_trans(proto1) - tri_trans([1,1,1]); 
    theta_proto1 = wrapTo2Pi(atan2(q(2), q(1))); 

    % find angle pointing from 1,1,1 to proto2
    q = tri_trans(proto2) - tri_trans([1,1,1]); 
    theta_proto2 = wrapTo2Pi(atan2(q(2), q(1))); 

    x = make_arch_in_triangle([1,1,1], proto1, n_cond, ...
        'angle_proto1', theta_proto1, ...
        'angle_proto2', theta_proto2, ...
        'n_cond_beyond_proto1', n_cond_beyond_proto1, ...
        'n_cond_beyond_proto2', n_cond_beyond_proto2, ...
        'direction', direction); 

    ratios = tri_trans_inv(x); 
    iois = ratios * pattern_period; 

    % project_ratios_triangle(x)
    % project_ratios_triangle(ratios)
    % project_ratios_triangle(iois)
    % 
    % 


else
    error('requested interpolation type %s not available...', interp_type); 
end










