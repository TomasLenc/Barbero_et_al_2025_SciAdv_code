function f = project_ratios_triangle(ratios, varargin)
% THis function projects a set of 3-interval rhythmic patterns onto a 2-D space
% (triangle). The input should be a [n_conditions x 3] matrix of interval
% ratios, where each row is a condition. 

parser = inputParser(); 

addParameter(parser, 'target_ratios', []); 

parse(parser, varargin{:}); 

target_ratios = parser.Results.target_ratios; 

%%

% get 3 corner points of the triangle
x1 = tri_trans([1 0 0]); 
x2 = tri_trans([0 1 0]); 
x3 = tri_trans([0 0 1]); 

f = figure('position',  1.5 * [675 311 785 651], 'color', 'white'); 

% plot the triangle 
h = fill([x1(1), x2(1), x3(1)], [x1(2), x2(2), x3(2)], [.5 .5 .5], 'FaceAlpha', 0.1); 
hold on 

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
    coords = tri_trans(int_ratios{i} / sum(int_ratios{i})); 
    scatter(coords(1), coords(2), 30, 'k', 'filled'); 
    text(coords(1), coords(2), ...
         sprintf(' %s', num2str(int_ratios{i}, '%d')),...
         'FontSize', 16); 
end

% plot user-requested target ratios
for i=1:length(target_ratios)
    coords = tri_trans(target_ratios{i} / sum(target_ratios{i})); 
    scatter(coords(1), coords(2), 30, 'bo', 'filled'); 
    text(coords(1), coords(2), ...
         sprintf(' %s', num2str(target_ratios{i}, '%d')),...
         'FontSize', 16, 'Color', 'b'); 
end

% plot the actual ratios for each condntion as a red point 
for i_cond=1:size(ratios, 1)
    if size(ratios, 2) == 3
        coords = tri_trans(ratios(i_cond, :) / sum(ratios(i_cond, :))); 
    else
        coords = ratios(i_cond, :); 
    end
    plot(coords(1), coords(2), 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r'); 
end


axis off




