function r = tri_trans_inv(x, pattern_duration)
% This function transforms from 2D coordinates on a triangle back to the 3D
% interval space of 3-interval rhythms
%
% Parameters
% ----------
% x : array_like, shape=[n_patterns,2]
%     X and Y coordinate of the pattern in the 2-D space where the triangle
%     lives. 
% pattern_duration : float, optional
%     Total duration of the output pattern. 
% 
% Returns 
% -------
% r :  array_like, shape=[n_patterns,3]
%     Vector with 3 elements corresponding to the duration of each interval
%     in the 3-interval pattern. 
%

x = ensure_row(x); 
r = nan(size(x,1), 3); 

for i=1:size(x, 1)

    r(i,3) = x(i,2) / sqrt(3/4);  
    r(i,2) = x(i,1) - 1/2*r(i,3); 
    r(i,1) = 1 - r(i,2) - r(i,3); 

    if nargin == 2
        r(i,:) = r(i,:) * pattern_duration; 
    end

end

end