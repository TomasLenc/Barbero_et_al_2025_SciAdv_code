function bounds_ratios = bound_to_ratio(ratio_per_condition, bounds)
% This fuction converts from boundary location in terms of index to an
% actual ratio value. 
%
% Parameters
% ----------
% ratio_per_condition :  array_like
%     Ratio value per condition. 
% bounds :  array_like
%     Boundary indices to convert
% 
% Returns 
% -------
% bounds_ratios : array_like
%     Boundary values converted to ratio values. 
% 

ratios_midpoints = ([ratio_per_condition(2:end); nan] + ratio_per_condition) / 2; 

bounds_ratios = ratios_midpoints(bounds); 


% bounds_ratios = ratio_per_condition(bounds+1); 
