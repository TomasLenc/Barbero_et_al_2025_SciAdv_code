function [f] = plot_best_model(rdm, model, r_obs, p, r_null_dist, varargin)
% Plot a summary of the best-fitting model. 
% 
% Parameters
% ----------
% rdm : array_like, shape=[n, n]
%     RDM of the data, comprising n conditions. 
% model : array_like, shape=[n, n]
%     RDM of the model. 
% r_obs : float
%     Observed correlation value between the rdm and model. 
% p : float
%     P-value of the correlation between the rdm and model. 
% r_null_dist : array of floats
%     Null distribution of the correlatino values obtained using
%     permutation testing. 
% z_snr : array of floats, optional
%     SNR values of the response across conditions. Mean and SD will be
%     plotted as text in teh figure. 
% 
% Returns 
% -------
% f : figure object 
%     Handle to the generated figure. 


parser = inputParser(); 

addParameter(parser, 'z_snr', []); 
addParameter(parser, 'flip', false); 
addParameter(parser, 'ratios', []); 
addParameter(parser, 'highlight_categs', []); 

parse(parser, varargin{:}); 

z_snr = parser.Results.z_snr; 
flip_matrix = parser.Results.flip; 
ratios = parser.Results.ratios; 
highlight_categs = parser.Results.highlight_categs; 

%% 

f = figure('pos', [680 801 402 161], 'color', 'white'); 

%%

ax = subplot(1, 3, 1); 

plot_rdm(rdm, 'ax', ax, 'flip', flip_matrix, ...
         'ratios', ratios, 'highlight_categs', highlight_categs); 

if ~isempty(z_snr)
    if length(z_snr) > 1
        tit = title(ax, sprintf('zSNR = %.2f \n+/- sd %.2f', mean(z_snr), std(z_snr))); 
    else
        tit = title(ax, sprintf('zSNR = %.2f', z_snr)); 
    end
end

%%

ax = subplot(1, 3, 2); 

plot_rdm(model, 'ax', ax, 'flip', flip_matrix,...
         'ratios', ratios, 'highlight_categs', highlight_categs); 
     
%%

ax = subplot(1, 3, 3); 

h = histogram(ax, r_null_dist, 'DisplayStyle', 'stairs', 'linew', 2);
hold on 
plot(ax, [r_obs, r_obs], [0, max(h.Values)], ':r', 'linew', 2)
ax.YAxis.Visible = 'off'; 
ax.XLim = [-1, 1]; 
axis square
box off

ax.Title.String = sprintf('r = %.2f (p = %.2g)', r_obs, p); 



