function f = plot_erp_overlay(erps, fs, varargin)
% Plot overlay of erps across conditions. Condition is indicated by color
% gradient. 
% 
% Parameters
% ----------
% erps : array_lik, shape [n_condiditons, time]
%     Array of time-domain erp data. Each row corresponds to one condition.
% fs : float
%     Sampling rate
% frex : array_like
%     Frequencies of interests. These will be highlighted in the plot. 
% linew :  float, optional, default=2
%     Line width. 
% ax :  axes object, optional
%     If passed, we will plot in these axes. If not, a new figure will be
%     opened and returned. 
% 
% Returns 
% -------
% f : figure object 
%     Generated figure. 
% 
parser = inputParser; 

addParameter(parser, 'ax', []); 
addParameter(parser, 'linew', 2); 
addParameter(parser, 'ratios', []); 

parse(parser, varargin{:}); 

ax = parser.Results.ax; 
linew = parser.Results.linew; 
ratios = parser.Results.ratios; 

if isempty(ax)
    f = figure('Color', 'white', 'Position', [680 794 329 184]); 
    ax = axes(f); 
else
    f = []; 
end

dur = size(erps, 2) / fs; 
t = [0 : size(erps, 2)-1] / fs; 
n_cond = size(erps, 1); 
% cols = brewermap(n_cond, 'PRGn'); 

cols = customcolormap([0, 1], {'#e87b0e', '#6a24a3'}, n_cond);


for i_cond=1:n_cond
  
    % plto the erp 
    plot(t, erps(i_cond, :), 'color', cols(i_cond, :), 'linew', linew)
    hold(ax, 'on'); 
            
end

box(ax, 'off'); 

cbar = colorbar; 
colormap(cols); 
cbar.Ticks = [0 : 2 : n_cond-1] / (n_cond-1); 
cbar.TickLabels = [1 : 2 : n_cond]; 

if ~isempty(ratios)
    ratio_labs = get_ratio_labs(ratios); 
    cbar.Ticks = [0, n_cond-1] / (n_cond-1); 
    cbar.TickLabels = ratio_labs([1, end]);     
end

ax.FontSize = 12; 
ax.XLim = [0, dur]; 
ax.XTick = [0, dur];

y_range = max(erps(:)) - min(erps(:)); 

ax.YLim = [min(erps(:)) - 0.1 * y_range, ...
           max(erps(:)) + 0.1 * y_range]; 
% ax.YTick = [-2, 0, 2]; 
