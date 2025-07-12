function [f, ax] = plot_rdm(rdm, varargin)
% Plot any RDM. 
% 
% Parameters
% ----------
% rdm : array_like, shape=[n, n]
%     Representational dissimiliary matrix with n dimensions. 
% ax :  axes object, optional
%     Handle to existing axes. If passed, we will plot the RDM there. 
% ratios : array_like, optional, shape=[1, n]
%     Ratio value for each condition that will be used as axis labels. 
% 
% Returns 
% -------
% f : figure object 
%     Handle to the generated figure. 
% ax : axes object 
%     Handle to the generated axes. 
% 
parser = inputParser(); 

addParameter(parser, 'ax', []); 
addParameter(parser, 'ratios', []); 
addParameter(parser, 'flip', false); 
addParameter(parser, 'highlight_categs', []); 

parse(parser, varargin{:}); 

ax = parser.Results.ax; 
ratios = parser.Results.ratios; 
flip_matrix = parser.Results.flip; 
categs = parser.Results.highlight_categs; 

if ~isempty(categs)
    assert(~isempty(ratios), ...
        'You asked to plot some categories as points but didnt provide interval ratios for each condition'); 
end

if flip_matrix
    rdm = rot90(rdm, 2)'; 
    ratios = flip(ratios); 
end

f = []; 
if isempty(ax)
    f = figure('color', 'white', 'Position',  [680 797 114 181]); 
    ax = axes(f); 
end

imagesc(ax, rdm);

% if there are any categories (i.e. prototype patterns) that should be
% shown as red points over the RDM, let's find their location 
categs_idx = nan(1, length(categs));
for i=1:length(categs)
    categ = categs{i} ./ sum(categs{i});
    [~, categs_idx(i)] = min(cellfun(@norm, num2cell(ratios - categ, 2))); 
end
hold(ax, 'on'); 
for i=1:length(categs_idx)
    plot(categs_idx(i), categs_idx(i), 'ro', 'MarkerFaceColor', 'r')
end

axis(ax, 'square'); 
ax.XAxisLocation = 'top'; 

% make sure the y-axis direction is reverse (the default differes depending
% on whether the axes are obtained from the builtin subplot function or
% from panel)
ax.YDir = 'reverse'; 

if ~isempty(ratios)
    ratio_labs = get_ratio_labs(ratios); 
    ax.XTick = [1, length(ratios)]; 
    ax.YTick = [1, length(ratios)]; 
    ax.XTickLabel = ratio_labs([1,end]); 
    ax.YTickLabel = ratio_labs([1,end]); 
else
    ax.XTick = []; 
    ax.YTick = []; 
    ax.XAxis.Visible = 'off'; 
    ax.YAxis.Visible = 'off'; 
end



