function [f, ax] = plot_r_bar(r_all, p_all, varargin)
% Plot a bar correponding to the mean correlation value, and individual
% correlation values as overlaid points. Plot significant points as filled,
% and non-significant as empty. 
% 
% Parameters
% ----------
% r_all : array_like, shape=[1, n]
%     Correlation values for n participants. 
% p_all : array_like, shape=[1, n]
%     P-values for n participants. 
% z_snr : array_like, shape=[1, n], optional
%     SNR values for n participants. If passed, this will be shown as
%     opacity of each individual point.
% color : array_like, shape=[1, 3], optional, default=[0 0.4470 0.7410]
%     RGB triplet for the color that will be used to plot.  
% y_maxlim : fload, optional, default=1
%     Maximum limit of the y axis. 
% 
% Returns 
% -------
% f : figure object 
%     Handle to the generated figure. 


parser = inputParser; 

addParameter(parser, 'ax', []); 
addParameter(parser, 'z_snr', []); 
addParameter(parser, 'y_maxlim', 1); 

parse(parser, varargin{:}); 

ax = parser.Results.ax; 
z_snr = parser.Results.z_snr; 
y_maxlim = parser.Results.y_maxlim; 

%%

c_group = [49, 97, 127] / 255; 
c_ind = [135, 175, 201]/255; 

f = []; 
if isempty(ax)
    f = figure('color', 'white', 'Position', [680 638 147 228]); 
    ax = axes(f); 
end

h = bar(ax, 0, mean(r_all), ...
        'FaceColor', 'none', ...
        'EdgeColor', c_group, ...
        'LineWidth', 3); 
    
hold(ax, 'on'); 

idx_signif = p_all < 0.05; 
    
h_signif = scatter(ax, 0.5 * (rand(1, sum(idx_signif))-0.5), r_all(idx_signif), ...
        50, c_ind, 'filled'); 
h_signif.MarkerFaceAlpha = 0.79; 
h_signif.MarkerEdgeAlpha = 0.79; 

h_nonsignif = scatter(ax, 0.5 * (rand(1, sum(~idx_signif))-0.5), r_all(~idx_signif), ...
        50, c_ind); 
h_nonsignif.MarkerFaceAlpha = 0.79; 
h_nonsignif.MarkerEdgeAlpha = 0.79; 

n_cols = 128; 
cmap = customcolormap([0,1], [c_group; [1, 1, 1]], n_cols); 
cbar = colorbar(ax); 
cbar.Location = 'eastoutside'; 
cbar.Visible = 'off'; 
cbar.Label.String = 'z snr'; 
cbar.Label.Position(1) = cbar.Label.Position(1) * 0.5; 

if ~isempty(z_snr)    
    
    z_snr_maxlim = norminv(1 - 1e-5); 
    edges = [-Inf, linspace(0, z_snr_maxlim, n_cols-1), Inf]; 
    c_idx = discretize(z_snr, edges); 

    h_signif.CData = cmap(c_idx(idx_signif), :); 
    h_nonsignif.CData = cmap(c_idx(~idx_signif), :); 
    
    colormap(cmap); 
    cbar.Limits = [0, 1]; 
    cbar.Ticks = [0, 1]; 
    cbar.TickLabels = [0, round(z_snr_maxlim, 2)]; 
    
    cbar.Visible = 'on'; 

end

sem = std(r_all) / sqrt(length(r_all)); 
ci = sem * norminv(1 - 0.025); 

plot(ax, [0,0], [mean(r_all)-ci, mean(r_all)+ci], ...
    'color', c_group, 'LineWidth', 4); 

box(ax, 'off'); 

ax.XLim = [-0.7, 0.7]; 
ax.XTick = []; 
ax.YLim = [0, y_maxlim]; 
ax.YTick = [0, y_maxlim]; 
ax.TickLength = [0,0];  



