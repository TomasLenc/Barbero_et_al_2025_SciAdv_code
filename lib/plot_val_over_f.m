function [f, ax] = plot_val_over_f(freq, val, varargin)

parser = inputParser; 

addParameter(parser, 'ax', []); 
addParameter(parser, 'p', []); 
addParameter(parser, 'color', [0, 0, 0] / 255); 

parse(parser, varargin{:}); 

ax = parser.Results.ax; 
c = parser.Results.color; 
p = parser.Results.p; 

f = []; 
if isempty(ax)
    f = figure('pos', [809 724 299 145], 'color', 'white'); 
    pnl = panel(f); 
    ax = pnl.select(); 
end

val_mu = mean(val, 1); 
val_sem = std(val, [], 1) / sqrt(size(val, 1)); 

h = fill([freq, flip(freq)], [val_mu - val_sem, flip(val_mu + val_sem)], c, ...
    'edgecolor', 'none', 'facealpha', 0.2); 
hold on 
plot(freq, val_mu, '-o', 'linew', 2, 'color', c, 'markerfacecolor', c, 'markersize', 5)

if ~isempty(p)
    mask_signif = p < 0.05; 
    plot(freq(mask_signif), val_mu(mask_signif) + val_sem(mask_signif) + 0.02, 'k*');  
end

% ax.XLim = [0, 60]; 
% ax.YLim = [0.3, 0.7]; 
% ax.YTick = [0.4, 0.6]; 
