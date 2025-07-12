% Summarize RSA results for ITI ratio responses. 

clear 

%%

experiment = '2ioi'; 

par = get_par('experiment', experiment); 

subjects = num2cell(par.subjects); 

feat_name_resp = 'ratio'; 
feat_name_stim = 'ratio';

response = 'tap-onsets'; 

run_stats = false; 

plot_big_figure = false; 

%%

subdirs = fullfile(...
                sprintf('response-%s', response)); 

save_path = fullfile(par.summary_path, subdirs); 
mkdir(save_path); 

fname_save_base = fullfile(save_path, ...
        sprintf('response-%s', ...
                         response)); 

%%

fname_tbl = sprintf('response-%s_feat-%s_featStim-%s.csv', ...
                         response, ...
                         feat_name_resp, ...
                         feat_name_stim); 

tbl = readtable(fullfile(par.rsa_path, subdirs, fname_tbl)); 

res = []; 

for i_sub=1:length(subjects)

    res(i_sub).sub = subjects{i_sub}; 
    sub_str = sub_num2str(res(i_sub).sub); 

    load_path_sub = fullfile(par.rsa_path, subdirs, sub_str); 

    fname = fullfile(load_path_sub, ...
        sprintf('%s_feat-%s_method-pearson_rdm.mat', sub_str, feat_name_resp)); 

    rdm = load(fname); 
    res(i_sub).rdm = rdm.rdm; 


    fname = fullfile(load_path_sub, ...
        sprintf('%s_feat-%s_featStim-%s_method-pearson_nCateg-2_perm.mat',...
                sub_str, feat_name_resp, feat_name_stim)); 

    perm_res = load(fname); 

    res(i_sub).best_model = perm_res.best_model; 
    res(i_sub).boundary = perm_res.all_models_info{1, 'bound1'}; 
    res(i_sub).r = perm_res.r_obs; 
    res(i_sub).p_perm = perm_res.p_perm; 

    [res(i_sub).tested_bounds, idx] = sort(perm_res.all_models_info.bound1); 
    res(i_sub).tested_bounds = ensure_row(res(i_sub).tested_bounds); 

    res(i_sub).r_all_models = ...
                        ensure_row(perm_res.all_models_info.r(idx)); 

    res(i_sub).r_null_dist = perm_res.r_null_dist; 

    mask = tbl.sub == res(i_sub).sub & tbl.n_categ_fit == 2; 

    % let's also load the mean tapped ratio for each subject 
    data = load(fullfile(par.deriv_path, 'preproc', sub_str, ...
             sprintf('%s_response-tap-ratios.mat', sub_str))); 

    res(i_sub).stim_ratios = data.stim_ratios; 
    res(i_sub).tap_ratios = data.tap_ratios; 
    res(i_sub).tap_ratios_all = data.tap_ratios_all; 
    res(i_sub).n_valid_tap_ratios = data.n_valid_tap_ratios; 

end

% save data 
save([fname_save_base, '_data.mat'], 'res'); 


tbl_ratios = cell2table(cell(0, 3), 'VariableNames', {'sub','stim_rat','tap_rat'}); 
for i_sub=1:length(res)
    for i_cond=1:par.n_cond
        new_row = [...
            {res(i_sub).sub}, ...
            {res(i_sub).stim_ratios(i_cond,1)}, ...
            {res(i_sub).tap_ratios(i_cond,1)}...
            ]; 
        tbl_ratios = [tbl_ratios; new_row]; 
    end
end
writetable(tbl_ratios, [fname_save_base, '_ratios.csv']); 


%% load stimulus 

load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 

stim_rdm = rsa_stim(par, [],'feat', 'ratio'); 



%% group-level permutation test

if exist([fname_save_base, '_permGroup.mat']) && ~run_stats
    clear res_group
    load([fname_save_base, '_permGroup.mat']); 
else
    res_group = []; 
    [res_group.r_obs, res_group.r_null_dist, res_group.p_perm] = ...
                                    perm_rdm_corr_group(...
                                              {res.rdm}, 2, 10000, ...
                                              'covariate', stim_rdm, ...
                                              'n_mini_categ_skip', 1); 

    save([fname_save_base, '_permGroup.mat'], 'res_group'); 
end       

%% grand average RDM 

load_path = fullfile(par.rsa_path, subdirs, 'sub-grand'); 
fname = sprintf('sub-grand_feat-%s_method-pearson_rdm', feat_name_resp); 
res_grand = load(fullfile(load_path, fname)); 


%% big summary figure 

if plot_big_figure

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


    % overlay of best fitting models across all subjects 
    plot_boundary_distr([res.boundary], par.n_cond, 'alpha', 0.1, 'pnl', pnl(2, 2, 1)); 
    pnl(2, 2, 1).title('all'); 

    % overlay of best fitting models across significant subjects 
    idx_signif = [res.p_perm] < 0.05; 
    plot_boundary_distr([res(idx_signif).boundary], par.n_cond, 'alpha', 0.1, 'pnl', pnl(2, 2, 2)); 
    pnl(2, 2, 2).title('signif categ'); 

    % distribution of r values across all tested models (i.e. across category
    % boundaries) 
    ax = pnl(2, 1, 2, 3).select(); 
    plot_r_over_all_models(cat(1, res.r_all_models), res(1).tested_bounds, 'ax', ax);

    % group-level permutation test 
    ax = pnl(2, 1, 1, 2).select(); 
    h = histogram(ax, res_group.r_null_dist, 'DisplayStyle', 'stairs', 'linew', 2);
    hold(ax, 'on'); 
    plot(ax, [res_group.r_obs, res_group.r_obs], [0, max(h.Values)], ':r', 'linew', 2)
    title(sprintf('r = %.2f, p = %.1g', res_group.r_obs, res_group.p_perm)); 
    xlabel('r'); 

    % grand average RDM 
    ax = pnl(2, 1, 1, 1).select(); 
    plot_rdm(res_grand.rdm, 'ax', ax); 
    pnl(2, 1, 1, 1).title('grand avg'); 

    % plot individual RDMs
    pnl_ind_rdm = pnl(1); 

    plot_all_rdms(pnl_ind_rdm, ...
        [res.sub], ...
        {res.rdm}, ...
        {res.best_model}, ...
        [res.p_perm] ...
        )

    pnl_ind_rdm.de.margin = [0, 0, 0, 8]; 
    pnl_ind_rdm.margin = [20, 2, 2, 15]; 

    pnl.marginleft = 15; 
    pnl.margintop = 15;
    pnl.fontsize = 12; 

    tit = pnl_ind_rdm.title('lf'); 
    tit.Position(2) = 1.03; 
    tit.FontWeight = 'bold';

    print(f, '-dpng', '-painters', '-r600', [fname_save_base, '_summary.png']); 

end





%% individual plots to assemble Figures for the paper
cmap_ratios = customcolormap([0, 1], {'#e87b0e', '#6a24a3'}, par.n_cond);


%% individual RDMs

f = figure('pos',[809 69 113 959], 'color', 'white'); 
pnl = panel(f); 
plot_all_rdms(pnl, ...
    [res.sub], ...
    {res.rdm}, ...
    {res.best_model}, ...
    [res.p_perm] ...
    )
pnl.margin = [10, 1, 1, 10]; 

print(f, [fname_save_base, '_indRDMs.eps'], '-depsc', '-painters')
close(f); 

%% FIGURE 2

% RHO of individual subjects 
f = figure('pos', [809 696 80 173], 'color', 'white'); 
pnl = panel(f); 
ax = pnl.select(); 

plot_r_bar([res.r], [res.p_perm], 'y_maxlim', 1, 'ax', ax); 

pnl.margin = [8,4,0,1]; 

saveas(f, [fname_save_base, '_rhoInd.svg']); 
close(f); 


% overlay of best fitting models across all significant subjects 
rdms = {res.best_model}; 
mask_signif = [res.p_perm] < 0.05; 

f = plot_rdms_overlay(rdms(mask_signif), 'alpha', 0.1); 

saveas(f, [fname_save_base, '_bestModelOverlaySignif.svg']); 
close(f); 


% median best-fitting model 
median_bound = median([res.boundary]); 
median_model = get_2categ_model(13, median_bound); 
f = figure('pos', [809 756 145 113], 'color', 'white'); 
pnl = panel(f); 
ax = pnl.select(); 
plot_rdm(median_model, 'ax', ax); 
saveas(f, [fname_save_base, '_medianBestModelRDM.svg']); 
close(f); 

rsa_bounds = [res.boundary]; 
median(rsa_bounds)
m = bootstrp(1000, @median, rsa_bounds);
prctile(m, [2.5, 100-2.5])

bounds_ratios = bound_to_ratio(stim.ratios(:, 1), [res.boundary]); 
median(bounds_ratios)
m = bootstrp(1000, @median, bounds_ratios);
prctile(m, [2.5, 100-2.5])



% grand avearge RDM
f = figure('pos', [809 756 145 113], 'color', 'white'); 
pnl = panel(f); 
ax = pnl.select(); 
plot_rdm(res_grand.rdm, 'ax', ax); 
pnl.margin = [1,1,1,1];

saveas(f, [fname_save_base, '_grandAvgRDM.svg']); 
close(f); 


% MDS            
cols = customcolormap([0, 1], {'#e87b0e', '#6a24a3'}, par.n_cond);
f = figure('pos', [606 688 165 181], 'color', 'white'); 

Y = {}; 
for i_sub=1:length(res)
    D = -res(i_sub).rdm + 1; 
    Y{i_sub} = cmdscale(D, 2);
%         if size(Y{i_sub}, 2) == 1
%             Y{i_sub}(:, 2) = 0; 
%         end
end

if all(cellfun(@(x) size(x, 2) == 2, Y))
    Y = mean(cat(3, Y{:}), 3); 
    s = scatter(Y(:, 1), Y(:, 2), 100, cols, 'filled'); 
    axis([min(Y(:,1)) - 0.2 * range(Y(:,1)), ...
          max(Y(:,1)) + 0.2 * range(Y(:,1)), ...
          min(Y(:,2)) - 0.2 * range(Y(:,2)), ...
          max(Y(:,2)) + 0.2 * range(Y(:,2))])
else
    Y = {}; 
    for i_sub=1:length(res)
        D = -res(i_sub).rdm + 1; 
        Y{i_sub} = cmdscale(D, 1);
    end
    Y = mean(cat(3, Y{:}), 3); 
    s = scatter(Y, repmat(0, 1, 13), 100, cols, 'filled'); 
    axis([min(Y(:,1)) - 0.2 * range(Y(:,1)), ...
          max(Y(:,1)) + 0.2 * range(Y(:,1)), ...
          -1, ...
          1])
end

axis square           
ax = gca; 
ax.XTick = []; 
ax.YTick = []; 
ax.XDir = 'reverse'; 

saveas(f, [fname_save_base, '_MDS.svg']); 
close(f); 

%% check number of valid data points 

min(cat(2, [res.n_valid_tap_ratios]), [], 1)
max(cat(2, [res.n_valid_tap_ratios]), [], 1)

mean(cat(2, [res.n_valid_tap_ratios]), 1)




%% sigmoid fit and ITI distribution 


f = figure('color', 'white', 'Position', [680 789 231 189]); 
pnl = panel(f); 
pnl.pack({80, 20}, {70, 30});

c_light = [135, 175, 201]/255; 
c_dark = [49, 97, 127]/255; 

%%%%%%%%%%%%%%%%%%  fit sigmoid to all data  %%%%%%%%%%%%%%%%%%

stim_ratios_all = repmat(stim.ratios(:, 1), length(res), 1); 
tap_ratios_all = cat(1, res.tap_ratios); 
tap_ratios_all = tap_ratios_all(:, 1); 

%     tap_ratios_all = cellfun(@(x) cat(1,x), {res.tap_ratios_all}, 'uni', 0); 
%     tap_ratios_all = cat(2, tap_ratios_all{:}); 
%     stim_ratios_all = repelem(stim.ratios(:, 1), sum(cellfun(@(x) size(x, 1), tap_ratios_all), 2)); 
%     tap_ratios_all = cellfun(@(x) cat(1,x{:}), num2cell(tap_ratios_all, 2), 'uni', 0); 
%     tap_ratios_all = cat(1, tap_ratios_all{:});         
%     tap_ratios_all = tap_ratios_all(:, 1); 

stim_ratios_all_for_pred = [...
                    min(stim_ratios_all) : ...
                    (max(stim_ratios_all)-min(stim_ratios_all))/100 : ...
                    max(stim_ratios_all)...
                    ];

[param, stats, y_pred] = fit_sigmoid(stim_ratios_all, tap_ratios_all, [], [], stim_ratios_all_for_pred); 

% "min", "max", "x50", "slope"
min_ = param(1); 
max_ = param(2); 
x50_ = param(3); 
slope_ = param(4); 

ax = pnl(1, 1).select(); 
hold(ax, 'on'); 

plot(stim_ratios_all, tap_ratios_all, '.', 'color', c_light)
h = plot(stim_ratios_all_for_pred, y_pred, '-', 'color', c_dark, 'linew', 3);

plot(ax, [1/2, 2/3], [1/2, 2/3], 'k')

x_lims = [min(stim_ratios_all)*0.99 max(stim_ratios_all)*1.01]; 
y_lims = [min(tap_ratios_all)*0.99 max(tap_ratios_all)*1.01]; 
xlim(x_lims)
ylim(y_lims)

ratio_labs = get_ratio_labs(stim.ratios); 

ax.XTick = []; 
yticks(stim.ratios([1, end], 1)); 
yticklabels(ratio_labs([1, end]))

%%%%%%%%%%%%%%%%%%  LOC linear vs. sigm model  %%%%%%%%%%%%%%%%%%

r2_lm = nan(1, length(res)); 
r2_sigm = nan(1, length(res)); 
for i_sub_lo=1:length(res)

    mask = ones(1, length(res), 'logical'); 
    mask(i_sub_lo) = 0; 

    x_train = repmat(stim.ratios(:, 1), length(res)-1, 1); 
    y_train = cat(1, res(mask).tap_ratios); 
    y_train = y_train(:, 1); 

    x_test = stim.ratios(:, 1); 

    y_obs = cat(1, res(~mask).tap_ratios); 
    y_obs = y_obs(:,1);         

    % linear model 
    param = polyfit(x_train, y_train, 1); 
    y_test = param(2) + param(1) * x_test; 
    rss = sum((y_obs - y_test).^2); 
    tss = sum(y_obs.^2); 
    r2_lm(i_sub_lo) = 1 - (rss/tss); 

    % sigmoid model 
    [param, stats, y_test] = fit_sigmoid(x_train, y_train, [], [], x_test); 
    rss = sum((y_obs - y_test).^2); 
    tss = sum(y_obs.^2); 
    r2_sigm(i_sub_lo) = 1 - (rss/tss); 

end

mean(r2_sigm)
mean(r2_lm)

ranksum(r2_sigm, r2_lm)

%%%%%%%%%%%%%%%%%%  bootstrap parameter values  %%%%%%%%%%%%%%%%%%

n_boot = 10000; 

min_boot_ = nan(1, n_boot); 
max_boot_ = nan(1, n_boot); 
x50_boot_ = nan(1, n_boot); 

x_boot = repmat(stim.ratios(:, 1), length(res), 1); 
for i_boot=1:n_boot
    if mod(i_boot, 100) == 0; fprintf('boot %d/%d\n', i_boot, n_boot); end
    % prepare boot data
    idx_sub_boot = randsample(length(res), length(res), true); 
    y_boot = cat(1, res(idx_sub_boot).tap_ratios); 
    y_boot = y_boot(:, 1); 

    % sigmoid model 
    param = fit_sigmoid(x_boot, y_boot, [], [], []); 
    min_boot_(i_boot) = param(1); 
    max_boot_(i_boot) = param(2); 
    x50_boot_(i_boot) = param(3); 
end

stim.ratios(:, 1) 

mean(x50_boot_)
prctile(x50_boot_, 2.5)
prctile(x50_boot_, 100-2.5)   


1 - mean(min_boot_ < 0.5)
mean(min_boot_)
prctile(min_boot_, 2.5)
prctile(min_boot_, 100-2.5)

1 - mean(max_boot_ < 2/3)
mean(max_boot_)
prctile(max_boot_, 2.5)
prctile(max_boot_, 100-2.5)

%%%%%%%%%%%%%%%%%% mean reproduction from edge conditions %%%%%%%%%%%%%%%%%%

tap_rat_1 = cellfun(@(x) x(1, 1) ,{res.tap_ratios}, 'uni', 1); 

bf10 = bf.ttest(tap_rat_1, 1/2)
[h, p, stat] = ttest(tap_rat_1, 1/2)

tap_rat_13 = cellfun(@(x) x(13, 1) ,{res.tap_ratios}, 'uni', 1); 
bf10 = bf.ttest(tap_rat_13, 2/3)
[h, p, stat] = ttest(tap_rat_13, 2/3)


%%%%%%%%%%%%%%%%%%  fit sigmoid to individual data  %%%%%%%%%%%%%%%%%%

min_ind_ = nan(1, length(res)); 
max_ind_ = nan(1, length(res)); 
x50_ind_ = nan(1, length(res)); 
for i_sub=1:length(res)
    stim_ratios_ind = stim.ratios(:, 1); 
    tap_ratios_ind = res(i_sub).tap_ratios(:, 1);

    stim_ratios_ind_for_pred = [min(stim_ratios_ind) : (max(stim_ratios_ind)-min(stim_ratios_ind))/100 : max(stim_ratios_ind)];
    [param, stats, y_pred] = fit_sigmoid(stim_ratios_ind, tap_ratios_ind, [], [], stim_ratios_ind_for_pred); 
    min_ind_(i_sub) = param(1); 
    max_ind_(i_sub) = param(2); 
    x50_ind_(i_sub) = param(3); 
end

bf10 = bf.ttest(min_ind_, 1/2)
[h, p, stat] = ttest(min_ind_, 1/2)

bf10 = bf.ttest(max_ind_, 2/3)
[h, p, stat] = ttest(max_ind_, 2/3)

ax = pnl(2, 1).select(); 
hold(ax, 'on'); 

% plot horizontal grey bars for each stimulus ratio 
for i_rat=1:size(stim.ratios, 1)    
    h = plot([stim.ratios(i_rat, 1), stim.ratios(i_rat, 1)], ...
         [-1, 2], 'k', 'linew', 3);
    h.Color(4) = 0.2; 
end

mu = mean(x50_ind_); 
sem = std(x50_ind_) / sqrt(length(res)); 
ci = sem * norminv(1 - 0.025); 

% plot individual estiamted thresholds 
plot([mu, mu], [-1, 2], 'color', c_dark, 'linew', 2); 
% plot([mu-ci, mu+ci], [0.5, 0.5], 'k', 'linew', 3); 

%     h = plot(x50_all_, rand(1, length(res)), 'o', 'MarkerEdgeColor', ...
%         'none', 'MarkerFaceColor', 'r', 'MarkerSize', 5); 
h = scatter(x50_ind_, rand(1, length(res)), 10, ...
            'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', c_light, ...
            'MarkerEdgeColor', 'none'); 

ax.XLim = x_lims; 
ax.YLim = [-1, 2]; 
ax.YTick = []; 
ax.XTick = stim.ratios([1, end], 1); 
ax.XTickLabel = ratio_labs([1, end]); 


%%%%%%%%%%%%%%%%%%  ITI distribution  %%%%%%%%%%%%%%%%%%

dt = 0.0001; 
fs = 1/dt; 

% hard-code sigma? 
sigma = 0.040; 

t = [0 : fs*1-1] / fs; 
t_rat = t ./ par.cycle_dur; 

tap_ratios_all = cat(1, res.tap_ratios); 
tap_ratios_all = tap_ratios_all(:,1); 

tap_iti_all = sort(tap_ratios_all) * par.cycle_dur; 
tap_int_idx = dsearchn(t', tap_iti_all); 
tap_int_vec = zeros(1, length(t)); 
tap_int_vec(tap_int_idx) = 1; 


% -------- Jacoby 2017 ---------

dens_ind = 1/(sqrt(2*pi)*sigma) * exp(-1/2 * ((tap_iti_all-t)/sigma).^2); 
dens = sum(dens_ind, 1) ./ length(tap_iti_all); 

%     figure
%     subplot 211
%     plot(t_rat, tap_int_vec)
%     xlim([1/2-0.2, 2/3+0.2])        
%     subplot 212
%     plot(t_rat, dens)
%     xlim([1/2-0.2, 2/3+0.2])

% this is a probability distribution!
sum(dens * dt)

% polot histogram 
n_bins = 30; % 30 when working with avergaed ITI, 90 with non-averaged?

ax = pnl(1, 2).select(); 
hold(ax, 'on'); 
h = histogram(ax, tap_ratios_all, n_bins); 
h.LineWidth = 2; 
h.FaceColor = c_light; 
h.EdgeColor = 'none'; 

ax.YAxis.Visible = 'off'; 
view(ax, [90 -90]) %// instead of normal view, which is view([0 90])
ax.XLim = y_lims; 
ax.XTick = stim.ratios([1, end], 1); 
ax.XTickLabel = []; 

% Chi square goodness of fit against unniform distribution  
counts = h.Values;     
N = sum(counts); 
expected_counts = repmat(N/n_bins, 1, n_bins); 
plot(h.BinEdges(1:end-1), expected_counts, '--', 'linew', 2)
chi_sq = sum((counts - expected_counts).^2 ./ expected_counts); 
df = n_bins - 1;     
p = chi2cdf(chi_sq, df, 'upper'); 
txt = text(ax, 2/3, 1, sprintf('?2(%d)=%.2f, p=%.1g', df, chi_sq, p)); 
txt.Rotation = -90;

% scale density for plotting 
dens_to_plot = dens ./ max(dens) * prctile(counts, 70); 
plot(ax, t_rat, dens_to_plot, 'color', c_dark, 'linew', 2)

pnl.de.margin = [10, 5, 5, 2]; 
pnl(1, 1).marginbottom = 2; 
pnl(2, 1).margintop = 0; 
pnl(1, 1).marginright = 3; 
pnl(2, 1).marginright = 3; 
pnl(1, 2).marginleft = 0; 
pnl(2, 2).marginleft = 0; 
pnl.margin = [10, 10 , 3, 3]; 


saveas(f, [fname_save_base, '_sigmFit.svg']); 



%% FIGURE 4 - FIXED MODELS

% group-level permutation test with fixed model 
if exist([fname_save_base, '_permGroupFixModel.mat']) && ~run_stats
    clear res_group
    load([fname_save_base, '_permGroupFixModel.mat']); 
else
    res_group = []; 
    [res_group.r_obs, res_group.r_null_dist, res_group.p_perm, res_group.models, res_group.models_info] = ...
                            perm_rdm_corr_group_fix_model(...
                                      {res.rdm}, 2,  par.n_perm_rsa_group, ...
                                      'covariate', stim_rdm, ...
                                      'n_mini_categ_skip', 1); 

    save([fname_save_base, '_permGroupFixModel.mat'], 'res_group'); 
end

% individual level permutation test with fixed model
if exist([fname_save_base, '_permIndFixModel.mat']) && ~run_stats
    clear res_ind_fix_model
    load([fname_save_base, '_permIndFixModel.mat']); 
else
    res_ind_fix_model = []; 
    for i_sub=1:length(subjects)
        i_sub
        [res_ind_fix_model(i_sub).r_obs, ...
         res_ind_fix_model(i_sub).r_null_dist, ...
         res_ind_fix_model(i_sub).p_perm, ...
         res_ind_fix_model(i_sub).models, ...
         res_ind_fix_model(i_sub).models_info, ...
         ] = ...
                            test_each_categ_model(...
                                      res(i_sub).rdm, 2,  par.n_perm_rsa_group, ...
                                      'covariate', stim_rdm, ...
                                      'n_mini_categ_skip', 1); 
    end
    save([fname_save_base, '_permIndFixModel.mat'], 'res_ind_fix_model'); 
end

f = figure('pos', [652 621 300 311], 'color', 'white'); 
pnl = panel(f); 
ax = pnl.select(); 

n_sub = length(res_ind_fix_model); 
n_models = length(res_group.r_obs); 

tmp = {res_ind_fix_model.r_obs}; 
r_ind = cat(2, tmp{:})';  

tmp = {res_ind_fix_model.p_perm}; 
p_ind = cat(1, tmp{:}); 

c_group = [49, 97, 127] / 255; 
c_ind = [135, 175, 201]/255; 

for i_model=1:n_models

    h = bar(ax, i_model, res_group.r_obs(i_model), ...
            'FaceColor', 'none', ...
            'EdgeColor', cmap_ratios(i_model, :), ...
            'LineWidth', 3); 

    hold(ax, 'on'); 

    idx_signif = p_ind(:, i_model) < 0.05; 

    h_signif = scatter(ax, i_model + 0.1 * (rand(1, sum(idx_signif))-0.5), ...
                       r_ind(idx_signif, i_model), ...
                       50, cmap_ratios(i_model, :), 'filled'); 

    h_nonsignif = scatter(ax, i_model + 0.1 * (rand(1, sum(~idx_signif))-0.5), ...
                          r_ind(~idx_signif, i_model), ...
                          50, cmap_ratios(i_model, :)); 

    h_signif.MarkerFaceAlpha = 0.6;                   
    h_nonsignif.MarkerFaceAlpha = 0.6;                   

    r_sem = std(r_ind(:, i_model)) / sqrt(n_sub); 
    r_ci = norminv(1 - 0.025) * r_sem; 

    plot([i_model, i_model], ...
         [res_group.r_obs(i_model) - r_ci, ...
          res_group.r_obs(i_model) + r_ci], ...
          'color', 'k', 'linew', 4)

    txt = text(i_model, max(r_ind(:,i_model))+0.1, sprintf('%.2g', res_group.p_perm(i_model))); 
    txt.HorizontalAlignment = 'center'; 

end

ax.YLim = [-1, 1]; 
ax.YTick = [-1, 0, 1]; 
ax.YAxis.TickDirection = 'out'; 
ax.XAxis.TickDirection = 'out'; 

ax.XTick = [1 : n_models]; 
ax.XTickLabel = sort(res_group.models_info.bound1); 


saveas(f, [fname_save_base, '_fixedModels.svg']); 

%%

 f = plot_r_over_all_models(cat(1, res.r_all_models), res(1).tested_bounds);

 saveas(f, [fname_save_base, '_fixedModelsDistr.svg']); 

 %%

f = figure('pos',[652 45 175 887], 'color', 'white'); 
pnl = panel(f); 
pnl.pack('v', length(res)); 
for i_sub=1:length(res)
    ax = pnl(i_sub).select(); 
    plot_r_over_all_models(cat(1, res(i_sub).r_all_models), res(1).tested_bounds, 'ax', ax);
    ax.XTickLabel = []; 
    ax.XAxis.Label.String = ''; 
    ax.YAxis.Label.String = ''; 
end
pnl.de.margin = [1, 1, 1, 0]; 

saveas(f, [fname_save_base, '_fixedModelsDistrInd.svg']); 



 %%
 close all

