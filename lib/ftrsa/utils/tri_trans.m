function x = tri_trans(r)
% Transofmration from 3D space of interval durations to a point on an
% equilateral triangle. To get the transform, just watch what happens with
% each basis in the 3-D space when it's transformed into the 2-D space
% where the triangle is... put it in a tranformation matrix M and you're
% done :)
%
% Parameters
% ----------
% r :  array_like, shape=[3,1]
%     Vector with 3 elements corresponding to the duration of each interval
%     in the 3-interval pattern. 
% 
% Returns 
% -------
% x : array_like, shape=[2,1]
%     X and Y coordinate of the pattern in the 2-D space where the triangle
%     lives. 
% 

r = ensure_row(r); 

assert(size(r,2) == 3, 'the input must be a 3D vector!')

M = [0, 0; 
     1, 0; 
     1/2, sqrt(3/4)]; 

 
x = nan(size(r,1), 2); 
for i=1:size(r,1) 
    x(i,:) = M' * (ensure_col(r(i,:)) ./ sum(r(i,:))); 
end

