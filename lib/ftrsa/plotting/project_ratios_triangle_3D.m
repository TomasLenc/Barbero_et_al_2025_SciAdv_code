function f = project_ratios_triangle_3D(ratios, varargin)
% This function projects a set of 3-interval rhythmic patterns onto a 3-D
% space. The input should be a [n_conditions x 3] matrix of interval
% ratios, where each row is a condition.

parser = inputParser(); 

addParameter(parser, 'target_ratios', []); 

parse(parser, varargin{:}); 

target_ratios = parser.Results.target_ratios; 

%%

% get 3 corner points of the triangle
x1 = ([1 0 0]); 
x2 = ([0 1 0]); 
x3 = ([0 0 1]); 

f = figure('position',  1.5 * [675 311 785 651], 'color', 'white'); 

% plot the triangle 
h = fill3([x1(1), x2(1), x3(1)], ...
         [x1(2), x2(2), x3(2)], ...
         [x1(3), x2(3), x3(3)], ...
         [.6, .6, .6], 'FaceAlpha', 1); 
     
hold on 

ax = gca; 
ax.View = [71.3000   23.6000]; 
ax.TickLength = [0, 0];
ax.FontSize = 12; 


% plot a set of small integer ratios to be able to find our way through the
% triangle 
int_ratios = {
    [0, 0, 1]
    [0, 1, 0]
    [1, 0, 0]
    [1, 1, 1]
    [3, 2, 3]
    [3, 3, 2]
    [2, 3, 3]
    [1, 1, 2]
    [2, 1, 1]
    [1, 2, 1]
    [3, 2, 2]
    [2, 3, 2]
    [2, 2, 3]
    [1, 1, 3]
    [1, 3, 1]
    [3, 1, 1]
    [1, 2, 2]
    [2, 1, 2]
    [2, 2, 1]
    [1, 3, 2]
    [1, 2, 3]
    [3, 1, 2]
    [3, 2, 1]
    [2, 3, 1]
    [2, 1, 3]
    };  

for i=1:length(int_ratios)
    coords = (int_ratios{i} / sum(int_ratios{i})); 
    scatter3(coords(1), coords(2), coords(3), 30, 'k', 'filled'); 
end

% plot user-requested target ratios
for i=1:length(target_ratios)
    coords = (target_ratios{i} / sum(target_ratios{i})); 
    scatter(coords(1), coords(2), coords(3), 30, 'ro', 'filled'); 
end


