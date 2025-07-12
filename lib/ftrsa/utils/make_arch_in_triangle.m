function x = make_arch_in_triangle(circle_center_3d, point_on_circle_3d, n_points, varargin)

parser = inputParser; 

addParameter(parser, 'angle_proto1', 0); 
addParameter(parser, 'angle_proto2', pi); 
addParameter(parser, 'n_cond_beyond_proto1', 0); 
addParameter(parser, 'n_cond_beyond_proto2', 0); 
addParameter(parser, 'direction', 'anticlockwise'); % clockwise or anticlockwise

parse(parser, varargin{:}); 

angle_proto1 = parser.Results.angle_proto1; 
angle_proto2 = parser.Results.angle_proto2; 
n_cond_beyond_proto1 = parser.Results.n_cond_beyond_proto1; 
n_cond_beyond_proto2 = parser.Results.n_cond_beyond_proto2; 
direction = parser.Results.direction; 

%%

circle_center = tri_trans(circle_center_3d); 
circle_radius = norm(tri_trans(point_on_circle_3d) - circle_center); 

x = nan(n_points, 2); 
c = 1; 

angle_from_proto1_to_2 = angle_proto2 - angle_proto1; 

if strcmp(direction, 'clockwise')
    if angle_from_proto1_to_2 > 0
        angle_from_proto1_to_2 = - mod(2 * pi - angle_from_proto1_to_2, 2*pi); 
    end
elseif strcmp(direction, 'anticlockwise')
    if angle_from_proto1_to_2 < 0
        angle_from_proto1_to_2 = mod(2 * pi + angle_from_proto1_to_2, 2*pi); 
    end
else 
    error('unknow direction type "%s"', direction); 
end

for i=(1-n_cond_beyond_proto1) : (n_points+n_cond_beyond_proto2) 
    
    theta = angle_proto1 + ...
            (i-1) * angle_from_proto1_to_2/(n_points-1); 
        
    x(c,1) = cos(theta) * circle_radius + circle_center(1); 
    x(c,2) = sin(theta) * circle_radius + circle_center(2); 
    c = c+1; 
    
end

% if the prototype 1 has larger angle than proto2, we need to save the
% conditions in the opposite order
if angle_proto1 > angle_proto2
    x = flip(x, 1); 
end




































