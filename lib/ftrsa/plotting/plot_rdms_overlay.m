function f = plot_rdms_overlay(rdms, varargin)
% Plot overlay of a set of RDMs. Intended use case is to plot an overlay of
% categorical models for a set of participants. 
% 
% Parameters
% ----------
% rdms : cell of array_like, shape=[n, n]
%     RDMs to be plotted in the overlay. 
% alpha : float, optional, detauls=0.1
%     Opacity of each RDM in the plot. From 1 which is no transparency at
%     all, to 0 which is complete transparency. 
% 
% Returns 
% -------
% f : figure object 
%     Handle to the generated figure. 

parser = inputParser; 

addParameter(parser, 'alpha', 0.1); 
addParameter(parser, 'ax', []); 

parse(parser, varargin{:}); 

opacity = parser.Results.alpha; 
ax = parser.Results.ax; 

%%

if isempty(ax)
    f = figure('color', 'white', 'Position',  [680 797 114 181]); 
    ax = axes(f); 
else
    f = []; 
end

for i=1:length(rdms)
    
    rdm = rdms{i}; 
    h = imagesc(ax, rdm); 
    h.AlphaData = opacity; 
    
    hold(ax, 'on'); 
    
end

cmap = brewermap(128, 'Greys'); 
colormap(ax, cmap); 

axis(ax, 'square'); 
ax.XTick = []; 
ax.YTick = []; 

% make sure the y-axis direction is reverse (the default differes depending
% on whether the axes are obtained from the builtin subplot function or
% from panel)
ax.YDir = 'reverse'; 

