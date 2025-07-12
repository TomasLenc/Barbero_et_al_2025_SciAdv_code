function bounds = ratio_to_bound(ratio_per_condition, ratios)
% This fuction converts from exact ratio value to an index of condition
% that is closest to that particular boundary.
%
% Parameters
% ----------
% ratio_per_condition :  array_like
%     Ratio value per condition. 
% ratios :  array_like
%     Ratio values to convert. 
% 
% Returns 
% -------
% bounds_ratios : array_like
%     Boundary values converted to ratio values. 
% 

ratios_midpoints = ([ratio_per_condition(2:end); nan] + ratio_per_condition) / 2; 

[~, bounds] = min(abs(ratios_midpoints - ratios)); 


