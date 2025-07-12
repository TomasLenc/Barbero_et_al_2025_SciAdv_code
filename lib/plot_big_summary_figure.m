function f = plot_big_summary_figure(par, res, res_group, rdm_grand, mX_data)

f = figure('pos', [272 1 984 1046], 'color', 'white'); 

pnl = panel(f); 

pnl.pack('h', [20, 50, 30]); 

pnl(2).pack('v', [25, 25, 50]); 

pnl(2, 1).pack('v', [50, 50]); 
pnl(2, 1, 1).pack('h', 3); 
pnl(2, 1, 2).pack('h', 3); 

pnl(2, 2).pack('v', 3); 

% pnl.select('all'); 

% bar plot with individual subjects overlaid as points 
ax = pnl(2, 1, 2, 1).select(); 
plot_r_bar([res.r], [res.p_perm], 'y_maxlim', 1, 'ax', ax); 
ax.YLabel.String = 'r'; 

% bar plot with zSNR of individual subjects overlaid as points 
ax = pnl(2, 1, 2, 2).select(); 
plot_r_bar([res.z_snr], [res.p_perm], 'y_maxlim', Inf, 'ax', ax); 
ax.YTick = [0, round(mean([res.z_snr]), 1), floor(max([res.z_snr])*10)/10]; 
ax.YLabel.String = 'z snr'; 


% overlay of best fitting models across all subjects 
plot_boundary_distr([res.boundary], par.n_cond, 'alpha', 0.1, 'pnl', pnl(2, 2, 1)); 
pnl(2, 2, 1).title('all'); 

% overlay of best fitting models across significant subjects 
idx_signif = [res.p_perm] < 0.05; 
plot_boundary_distr([res(idx_signif).boundary], par.n_cond, 'alpha', 0.1, 'pnl', pnl(2, 2, 2)); 
pnl(2, 2, 2).title('signif categ'); 

% overlay of best fitting models across subjects with significant z-snr
idx_signif = [res.z_snr] > norminv(1 - 1e-3); 
plot_boundary_distr([res(idx_signif).boundary], par.n_cond, 'alpha', 0.1, 'pnl', pnl(2, 2, 3)); 
pnl(2, 2, 3).title('signif zsnr'); 
pnl(2, 2, 3, 2).xlabel('boundary location'); 


% distribution of r values across all tested models (i.e. across category
% boundaries) 
ax = pnl(2, 1, 2, 3).select(); 
plot_r_over_all_models(cat(1, res.r_all_models), res(1).tested_bounds, 'ax', ax);

% relationhsip between z-snr and r
ax = pnl(2, 1, 1, 3).select(); 
plot(ax, [res.z_snr], [res.r], 'ko')
axis(ax, 'square'); 
box(ax, 'off'); 
ax.XLabel.String = 'z snr'; 
ax.YLabel.String = 'r'; 

% group-level permutation test 
ax = pnl(2, 1, 1, 2).select(); 
h = histogram(ax, res_group.r_null_dist, 'DisplayStyle', 'stairs', 'linew', 2);
hold(ax, 'on'); 
plot(ax, [res_group.r_obs, res_group.r_obs], [0, max(h.Values)], ':r', 'linew', 2)
title(sprintf('r = %.2f, p = %.1g', res_group.r_obs, res_group.p_perm)); 
xlabel('r'); 

% grand average RDM 
ax = pnl(2, 1, 1, 1).select(); 
plot_rdm(rdm_grand, 'ax', ax); 
pnl(2, 1, 1, 1).title('grand avg'); 

% grand average mX 
pnl_mX = pnl(3); 
plot_mX_lw(mX_data.header_mX_to_plot, mX_data.data_mX_to_plot, ...
               [1 : mX_data.header_mX_to_plot.datasize(2)], ...
               mX_data.frex, ...
               'z_snr', mX_data.z_snr, ...
               'fmax', max([mX_data.frex, 15]), ...
               'pnl', pnl_mX); 

% plot individual RDMs
pnl_ind_rdm = pnl(1); 

plot_all_rdms(pnl_ind_rdm, ...
    [res.sub], ...
    {res.rdm}, ...
    {res.best_model}, ...
    [res.p_perm], ...
    [res.z_snr] ...
    )

pnl_ind_rdm.de.margin = [0, 0, 0, 8]; 
pnl_ind_rdm.margin = [20, 2, 2, 15]; 

pnl.marginleft = 15; 
pnl.margintop = 15;
pnl.fontsize = 12; 

tit = pnl_ind_rdm.title('lf'); 
tit.Position(2) = 1.03; 
tit.FontWeight = 'bold';




