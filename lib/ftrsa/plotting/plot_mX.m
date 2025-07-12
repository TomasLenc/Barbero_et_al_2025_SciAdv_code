function [f, pnl] = plot_mX(mX, freq, frex, varargin)
% Plot magnitude spectrum for each condition.  
% 
% Parameters
% ----------
% mX : cell array
%     Ech element contains the magnitude spectrum for one condition. 
% freq : array_like
%     Vector of frequencies. Must have the same size as each element in mX.
% frex : array_like
%     Frequencies of interests. These will be highlighted in the plot. 
% z_snr :  array_like, optional
%     SNR value for each condition. If passed, it will be plotted as text
%     label above each spectrum. 
% fmin :  float, optional, default=0.3
%     Minimum frequency limit of the plot. 
% fmax :  float, optional, default=max(frex)+1
%     Maximum frequency to plot. 
% 
% Returns 
% -------
% f : figure object 
%     Generated figure. 
% pnl : panel object 
%     Panel object in the figure. 
% 
parser = inputParser(); 

addParameter(parser, 'pnl', []); 
addParameter(parser, 'z_snr', []); 
addParameter(parser, 'fmin', 0.3); 
addParameter(parser, 'fmax', max(frex) + 1); 
addParameter(parser, 'show_x_line', false); 

parse(parser, varargin{:}); 

pnl = parser.Results.pnl; 
z_snr = parser.Results.z_snr; 
fmin = parser.Results.fmin; 
fmax = parser.Results.fmax; 
show_x_line = parser.Results.show_x_line; 

%%

n_cond = length(mX);

f = []; 
if isempty(pnl)
    f = figure('color', 'white', 'position', [631 71 355 70*n_cond]); 
    pnl = panel(f); 
end
pnl.pack('v', n_cond); 

pnl.de.margin = 0;

y_min = Inf; 
y_max = -Inf; 

frex_idx = dsearchn(ensure_col(freq), ensure_col(frex))'; 

min_freq_idx = dsearchn(ensure_col(freq), fmin); 
max_freq_idx = dsearchn(ensure_col(freq), fmax); 

for i_cond=1:n_cond              
        
    mX_to_plot = ensure_row(squeeze(mX{i_cond})); 
    
    ax = pnl(i_cond).select(); 

    plot(freq(min_freq_idx:max_freq_idx), ...
         mX_to_plot(min_freq_idx:max_freq_idx), ...
         'color', [49, 97, 127] / 255, ...
         'linew', 1); 
     
    hold on 
        
    plot(freq(frex_idx), mX_to_plot(frex_idx), 'ro'); 
    
    ax.XLim = [fmin, fmax]; 
    
    if ~isempty(z_snr)
        tit = title(ax, sprintf('zSNR = %.2f', z_snr(i_cond))); 
    end
    
    if i_cond ~= n_cond
        ax.XTick = []; 
        ax.YTick = []; 
        if ~show_x_line
            ax.XAxis.Visible = 'off'; 
        end
    end
        
    y_min = min(y_min, min(mX_to_plot(min_freq_idx:max_freq_idx))); 
    y_max = max(y_max, max(mX_to_plot(min_freq_idx:max_freq_idx))); 
    
end

for i_cond=1:n_cond
    ax = pnl(i_cond).select(); 
    ax.YLim = [y_min, y_max]; 
    ax.Title.Position(2) = ax.Title.Position(2) - ax.Title.Position(2) * 0.3; 
end

ax.YTick = [y_min, y_max]; 
ax.YTickLabel = unique(sort([0, round(y_min * 100)/100, round(y_max * 100)/100])); 
ax.TickDir = 'out'; 

pnl.de.margintop = 100/n_cond;

pnl.marginbottom = 10; 

pnl.margintop = 10; 

pnl.fontsize = 12; 


