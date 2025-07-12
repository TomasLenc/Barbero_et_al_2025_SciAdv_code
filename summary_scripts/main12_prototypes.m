% Perform the prototype analysis 

clear 


experiment = '2ioi'; 
par = get_par('experiment', experiment); 
subjects = par.subjects; 

par_proto = []; 

par_proto.response = 'eeg'; 
par_proto.harmonics = [1 : 6];  
par_proto.bound = 5; 

% par_proto.response = 'tap-force'; 
% par_proto.harmonics = [1 : 12];  
% par_proto.bound = 4; 

% par_proto.response = 'tap-impulse'; 
% par_proto.harmonics = [1 : 12];  
% par_proto.bound = 4; 

par_proto.feat_name = 'mX'; 

do_save = true; 

redo_boot = true; 

%%

par_proto.n_steps_between_proto = 37; 

par_proto.n_steps_beyond_proto1 = 0; 
par_proto.n_steps_beyond_proto2 = 39; % 10 takes us to 0.8 ratio (with granularity 13), 20 with granularity 26

load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 

par_proto.proto_ir = 1; 
par_proto.proto_ir_name = 'dirac'; 

[s_proto, ratios_proto] = make_slice_signals(stim.fs, ...
                                 par_proto.proto_ir, ...
                                 30, par.cycle_dur, ...
                                 'n_cond_beyond_proto1', par_proto.n_steps_beyond_proto1, ...
                                 'n_cond_beyond_proto2', par_proto.n_steps_beyond_proto2, ...
                                 'n_cond', par_proto.n_steps_between_proto, ...
                                 'proto1', [1,1],...
                                 'proto2', [2,1]); 
                             
                             
stim.ratios = stim.ratios(:, 1); 

ratios_proto = ratios_proto(:,1); 

ratios_proto_labs = get_ratio_labs(ratios_proto); 

stim_rat_step = stim.ratios(2) - stim.ratios(1); 

proto_ratio_step = ratios_proto(2) - ratios_proto(1); 

max_idx_hist_edges = [...
    ratios_proto(1)-proto_ratio_step/2 : ...
    proto_ratio_step : ...
    ratios_proto(end)+proto_ratio_step/2 ...
    ]; 

n_prototypes = length(s_proto); 

% find which prototypes correspond to stimulus ratios (note that these
% won't be exact)
cond_target_proto_idx = dsearchn(ratios_proto, stim.ratios)'; 

% get paths
save_path = fullfile(par.deriv_path, 'prototypes'); 
if ~isdir(save_path)
    mkdir(save_path)
end

%% analyze protoytpes

frex = 1/par.cycle_dur * par_proto.harmonics; 
harm_str = sprintf('harm-%s', vec_to_str(par_proto.harmonics, 'sep', '-')); 

fname_base = sprintf('response-%s_feat-%s_roi-front_%s_protoIR-%s_protoRes-%d_nProtoBeyond1-%d_nProtoBeyond2-%d', ...
                     par_proto.response, par_proto.feat_name, harm_str, par_proto.proto_ir_name, par_proto.n_steps_between_proto, ...
                     par_proto.n_steps_between_proto, par_proto.n_steps_beyond_proto2); 
                 
% get the protytopes for this seletion of par_proto.harmonics   
[X, freq] = get_X(s_proto, stim.fs); 
feat_proto = get_feat_from_X(X, freq, frex, 'feat', par_proto.feat_name); 

% try to find which prototype has zero magnitude at least in one harmonic                                     
[amp, idx] = sort(min(get_feat_from_X(X, freq, frex, 'feat', 'mX'), [], 1)); 
rat_proto_near_zero = unique(ratios_proto(idx(1:15)))'; 
idx_proto_with_zero = dsearchn(ratios_proto, rat_proto_near_zero'); 


%% plot proto info 

if n_prototypes < 40
    [f,pnl] = plot_time_mX_sim(s_proto, stim.fs, frex, 'filter_erp', false,...
                         'cond_labs', ratios_proto_labs); 
    f.Position =  [452 1208 610 1290]; 
    pnl.de.margintop = 6; 
    if do_save
        print(fullfile(save_path, [fname_base, '_proto-time-mX.png']), ...
                       '-dpng', '-painters', '-r600', f); 
        print(fullfile(save_path, [fname_base, '_proto-time-mX.eps']), ...
                       '-depsc', '-painters', '-r600', f); 
    end
end

f = figure('color', 'white', 'position', [631 71 600 930]); 
pnl = panel(f); 
ax = pnl.select(); 
pnl.de.margin = 0;
tmax = 5; 
N = round(tmax*stim.fs); 
t = [0:N-1]/stim.fs; 
for i_proto=1:n_prototypes
    plot(t, s_proto{i_proto}(1:N) + (n_prototypes - 1.5*i_proto), 'k', 'linew', 2); 
    hold(ax, 'on'); 
end
ax.YTick = flip(n_prototypes - 1.5*[1:n_prototypes]); 
ax.YTickLabel = flip(ratios_proto_labs); 
if do_save
    print(fullfile(save_path, [fname_base, '_proto-time.eps']), ...
          '-depsc', '-painters', '-r600', f); 
end
close(f);


%% get mX for each participant
    
X_all = cell(1, length(subjects)); 
mX_all = cell(1, length(subjects)); 
for i_sub=1:length(subjects)
    sub = subjects(i_sub); 
    sub_str = sub_num2str(sub);   
    % load data 
    fpath = fullfile(par.deriv_path, 'fft', sub_str); 
    fname = sprintf('sub-%03d_response-%s_concatEp-false_snr-none_X', sub, par_proto.response); 
    [header, data] = CLW_load(fullfile(fpath, fname));
    if strcmp(par_proto.response, 'eeg')
        % average frontocentral 
        chan_idx = get_chan_idx(header, par.roi_front_chans); 
        data = mean(data(:, chan_idx, 1, 1, 1, :), 2); 
    end
    % get X
    X_all{i_sub} = data; 
    % get noise subtracted mX
    freq = [0 : header.datasize(end)-1] * header.xstep; 
    mX_all{i_sub} = subtract_noise_bins(abs(data), ...
                                        par.snr_bins(1), par.snr_bins(2)); 
end


%% find distribution of winning prototype from bootstrapped grand average

par_proto.n_boot_level1 = 1000; 

max_idx_boot_obs = nan(par_proto.n_boot_level1, par.n_cond); 

parfor i_boot=1:par_proto.n_boot_level1

    if mod(i_boot, 10) == 0
        fprintf('%d/%d\n', i_boot, par_proto.n_boot_level1); 
    end
    
    idx_sub_boot = randsample(length(mX_all), length(mX_all), true); 

    if strcmp(par_proto.feat_name, 'X')
        X = squeeze(mean(cat(7, X_all{idx_sub_boot}), 7)); 
        feat_resp = get_feat_from_X(X, freq, frex, 'feat', 'X'); 
    elseif strcmp(par_proto.feat_name, 'mX')
        mX = squeeze(mean(cat(7, mX_all{idx_sub_boot}), 7)); 
        feat_resp = get_feat_from_X(mX, freq, frex, 'feat', 'mX'); 
    end

    [~, ~, ~, max_idx_boot_obs(i_boot,:)] = get_pariwise_distances(feat_proto, feat_resp, ...
                                         'normalize_columns', false, ...
                                         'plot', false); 
end

% find peaks on pre-par_proto.bound histogram 
rat_data = ratios_proto(reshape(max_idx_boot_obs(:, 1:par_proto.bound), [], 1)); 
x_kde = [ratios_proto(1) : 0.001 : ratios_proto(end)]; 
kde = ksdensity(rat_data, x_kde, 'bandwidth', 0.005); 
[pks, par_proto.peak_locs_pre] = findpeaks(kde, x_kde, 'MinPeakHeight', max(kde)*0.10); 


% find peaks on post-par_proto.bound histogram 
rat_data = ratios_proto(reshape(max_idx_boot_obs(:, par_proto.bound+1:end), [], 1)); 
x_kde = [ratios_proto(1) : 0.001 : ratios_proto(end)]; 
kde = ksdensity(rat_data, x_kde, 'bandwidth', 0.005); 
[pks, par_proto.peak_locs_post] = findpeaks(kde, x_kde, 'MinPeakHeight', max(kde)*0.10); 

par_proto.peak_locs_pre
par_proto.peak_locs_post


%% bootstrap the bootrap summaries

par_proto.n_boot_level2 = 300; 


if exist(fullfile(save_path, [fname_base, '_maxIdxGrandBoot_kde.mat'])) && ...
   ~redo_boot
    
load(fullfile(save_path, [fname_base, '_maxIdxGrandBoot_kde.mat']))

ratios_proto = ratios_proto(:,1); 
stim.ratios = stim.ratios(:,1); 
    
else

% on-off index calcualted for each condition by cetering the "on" window on
% the current condition ratio and taking all other ratios outside that
% window as "off". 
on_off_index_stim_ratio_boot = nan(par_proto.n_boot_level2, par.n_cond); 

% For the pre and post boundary histogram, we take a window approach to
% test whether the peaks are significant. For each peak, we will take a
% local surrounding window, and see if the density is higher inside vs.
% directly outside of the window. 
par_proto.ratios_to_test_pre = ratios_proto(dsearchn(ratios_proto, par_proto.peak_locs_pre')); 
par_proto.ratios_to_test_post = ratios_proto(dsearchn(ratios_proto, par_proto.peak_locs_post')); 

on_off_index_pre_bound_ratio_boot = nan(par_proto.n_boot_level2, length(par_proto.ratios_to_test_pre)); 
on_off_index_post_bound_ratio_boot = nan(par_proto.n_boot_level2, length(par_proto.ratios_to_test_post)); 

% get the x-axis limits of each on-peak window  
par_proto.on_lims_pre = nan(length(par_proto.ratios_to_test_pre), 2); 
par_proto.on_lims_post = nan(length(par_proto.ratios_to_test_post), 2);

% get the x-axis limits of each off-peak window  
par_proto.off_lims_below_pre = nan(length(par_proto.ratios_to_test_pre), 2); 
par_proto.off_lims_above_pre = nan(length(par_proto.ratios_to_test_pre), 2); 
par_proto.off_lims_below_post = nan(length(par_proto.ratios_to_test_post), 2); 
par_proto.off_lims_above_post = nan(length(par_proto.ratios_to_test_post), 2); 

for i_test_rat=1:length(par_proto.ratios_to_test_pre)
    par_proto.on_lims_pre(i_test_rat,1) = par_proto.ratios_to_test_pre(i_test_rat) - stim_rat_step/2; 
    par_proto.on_lims_pre(i_test_rat,2) = par_proto.ratios_to_test_pre(i_test_rat) + stim_rat_step/2; 

    par_proto.off_lims_below_pre(i_test_rat,1) = par_proto.ratios_to_test_pre(i_test_rat) - stim_rat_step; 
    par_proto.off_lims_below_pre(i_test_rat,2) = par_proto.ratios_to_test_pre(i_test_rat) - stim_rat_step + stim_rat_step/2; 

    par_proto.off_lims_above_pre(i_test_rat,1) = par_proto.ratios_to_test_pre(i_test_rat) + stim_rat_step - stim_rat_step/2; 
    par_proto.off_lims_above_pre(i_test_rat,2) = par_proto.ratios_to_test_pre(i_test_rat) + stim_rat_step; 
end

for i_test_rat=1:length(par_proto.ratios_to_test_post)
    par_proto.on_lims_post(i_test_rat,1) = par_proto.ratios_to_test_post(i_test_rat) - stim_rat_step/2; 
    par_proto.on_lims_post(i_test_rat,2) = par_proto.ratios_to_test_post(i_test_rat) + stim_rat_step/2; 

    par_proto.off_lims_below_post(i_test_rat,1) = par_proto.ratios_to_test_post(i_test_rat) - stim_rat_step; 
    par_proto.off_lims_below_post(i_test_rat,2) = par_proto.ratios_to_test_post(i_test_rat) - stim_rat_step + stim_rat_step/2; 

    par_proto.off_lims_above_post(i_test_rat,1) = par_proto.ratios_to_test_post(i_test_rat) + stim_rat_step - stim_rat_step/2; 
    par_proto.off_lims_above_post(i_test_rat,2) = par_proto.ratios_to_test_post(i_test_rat) + stim_rat_step; 
end

for i_boot_outer=1:par_proto.n_boot_level2

    % bootstrap mean of winning prototype determined from grand-avearge
    % -----------------------------------------------------------------
    
    if mod(i_boot_outer, 2) == 0
        fprintf('\n\nboot outer %d/%d\n', i_boot_outer, par_proto.n_boot_level2); 
    end        
    
    max_idx_boot_boot = nan(par_proto.n_boot_level1, par.n_cond); 
    parfor i_boot=1:par_proto.n_boot_level1
        if mod(i_boot, 100) == 0
            fprintf('%d/%d\n', i_boot, par_proto.n_boot_level1); 
        end
        idx_sub_boot = randsample(length(mX_all), length(mX_all), true); 
        if strcmp(par_proto.feat_name, 'X')
            X = squeeze(mean(cat(7, X_all{idx_sub_boot}), 7)); 
            feat_resp = get_feat_from_X(X, freq, frex, 'feat', 'X'); 
        elseif strcmp(par_proto.feat_name, 'mX')
            mX = squeeze(mean(cat(7, mX_all{idx_sub_boot}), 7)); 
            feat_resp = get_feat_from_X(mX, freq, frex, 'feat', 'mX'); 
        end
        [~, ~, ~, max_idx_boot_boot(i_boot,:)] = get_pariwise_distances(feat_proto, feat_resp, ...
                                             'normalize_columns', false, ...
                                             'plot', false); 
    end
    

    % find ON-OFF difference for window around stimulus ratio in ths condition
    % ----------------------------------------------------------------------

    for i_cond=1:par.n_cond
        
        if i_cond==1
            on_lim_l = stim.ratios(i_cond); 
            on_lim_u = stim.ratios(i_cond+1) - stim_rat_step/2; 
        elseif i_cond==par.n_cond
            on_lim_l = stim.ratios(i_cond-1) + stim_rat_step/2; 
            on_lim_u = stim.ratios(i_cond); 
        else                
            on_lim_l = stim.ratios(i_cond-1) + stim_rat_step/2; 
            on_lim_u = stim.ratios(i_cond+1) - stim_rat_step/2; 
        end

        rat_data = ratios_proto(max_idx_boot_boot(:, i_cond)); 

        n_on = sum((rat_data > on_lim_l) & (rat_data < on_lim_u));
        p_on = n_on / sum((ratios_proto > on_lim_l) & (ratios_proto < on_lim_u)); 

        n_off = sum((rat_data < on_lim_l) | (rat_data > on_lim_u));
        p_off = n_off / sum((ratios_proto < on_lim_l) | (ratios_proto > on_lim_u));

        on_off_index_stim_ratio_boot(i_boot_outer, i_cond) = p_on - p_off; 

    end    
    
    % pre boundary 
    % ------------
       
    rat_data = ratios_proto(reshape(max_idx_boot_boot(:, 1:par_proto.bound), [], 1)); 

    for i_test_rat=1:length(par_proto.ratios_to_test_pre)

        n_on = sum((rat_data > par_proto.on_lims_pre(i_test_rat,1)) & ...
                   (rat_data < par_proto.on_lims_pre(i_test_rat,2)));

        p_on = n_on / sum((ratios_proto > par_proto.on_lims_pre(i_test_rat,1)) & ...
                          (ratios_proto < par_proto.on_lims_pre(i_test_rat,2))); 

        n_off = sum(((rat_data > par_proto.off_lims_below_pre(i_test_rat,1)) & (rat_data < par_proto.off_lims_below_pre(i_test_rat,2))) | ...
                    ((rat_data > par_proto.off_lims_above_pre(i_test_rat,1)) & (rat_data < par_proto.off_lims_above_pre(i_test_rat,2))) );

        p_off = n_off / sum(((ratios_proto > par_proto.off_lims_below_pre(i_test_rat,1)) & (ratios_proto < par_proto.off_lims_below_pre(i_test_rat,2))) |...
                            ((ratios_proto > par_proto.off_lims_above_pre(i_test_rat,1)) & (ratios_proto < par_proto.off_lims_above_pre(i_test_rat,2)))...
                            ); 

        on_off_index_pre_bound_ratio_boot(i_boot_outer, i_test_rat) = p_on - p_off; 
        
    end    

    % post boundary 
    % -------------
    
    rat_data = ratios_proto(reshape(max_idx_boot_boot(:, par_proto.bound+1:end), [], 1), 1); 

    for i_test_rat=1:length(par_proto.ratios_to_test_post)

        n_on = sum((rat_data > par_proto.on_lims_post(i_test_rat,1)) & ...
                   (rat_data < par_proto.on_lims_post(i_test_rat,2)));

        p_on = n_on / sum((ratios_proto > par_proto.on_lims_post(i_test_rat,1)) & ...
                          (ratios_proto < par_proto.on_lims_post(i_test_rat,2))); 

        n_off = sum(((rat_data > par_proto.off_lims_below_post(i_test_rat,1)) & (rat_data < par_proto.off_lims_below_post(i_test_rat,2))) | ...
                    ((rat_data > par_proto.off_lims_above_post(i_test_rat,1)) & (rat_data < par_proto.off_lims_above_post(i_test_rat,2))) );

        p_off = n_off / sum(((ratios_proto > par_proto.off_lims_below_post(i_test_rat,1)) & (ratios_proto < par_proto.off_lims_below_post(i_test_rat,2))) |...
                            ((ratios_proto > par_proto.off_lims_above_post(i_test_rat,1)) & (ratios_proto < par_proto.off_lims_above_post(i_test_rat,2)))...
                            ); 

        on_off_index_post_bound_ratio_boot(i_boot_outer, i_test_rat) = p_on - p_off; 
        
    end    
    
end

end


%% plot 

show_stim_ticks = true; 
highlight_ratios = [0.5, 2/3, 3/4]; 

f = figure('color', 'w', 'pos', [571 365 270 452]); 

pnl = panel(f); 
pnl.pack('v', [15, 70, 15]); 
pnl(2).pack('v', size(max_idx_boot_obs, 2)); 
pnl.de.margin = 1; 


x_kde = [0 : 0.001 : 1]; 

lighten_color = @(col, factor) col + (1 - col) * factor;

% convert from index to ratio 
max_ratio_boot = ratios_proto(max_idx_boot_obs); 

proto_ratio_step = ratios_proto(2) - ratios_proto(1); 

hist_edges = [...
    ratios_proto(1)-proto_ratio_step/2 : ...
    proto_ratio_step : ...
    ratios_proto(end)+proto_ratio_step/2 ...
    ]; 

for i_cond=1:13
    
    if i_cond<=par_proto.bound
        col_light = [135, 175, 201]/255; 
        col_dark = [49, 97, 127] / 255; 
    else
        col_light = [191, 55, 21]/255; 
        col_dark = [168, 11, 0]/255; 
    end
    
    % get kde
    kde = ksdensity(max_ratio_boot(:,i_cond), x_kde, 'bandwidth', 0.005); 
    
    ax = pnl(2, i_cond).select();
    hold(ax, 'on'); 
    
    if show_stim_ticks       
         plot([stim.ratios(i_cond), stim.ratios(i_cond)], ...
             [0, max(kde)], ...
             '-', 'linew', 5, 'color', [112, 112, 112]/255); 
    end
    
    fill([x_kde, fliplr(x_kde)], [kde, zeros(size(kde))],  lighten_color(col_dark, 0.5), ...
         'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(x_kde, kde, 'color', lighten_color(col_dark, 0.5), 'linew', 2)

    ax.YAxis.Visible = 'off'; 
    ax.XLim = [hist_edges(1), hist_edges(end)];
    ax.XTick = []; 
end

% marginal for caterogy 1
col_light = [135, 175, 201]/255; 
col_dark = [49, 97, 127] / 255; 

ax = pnl(1).select(); 
data_for_hist = reshape(max_ratio_boot(:, 1:par_proto.bound), [], 1); 
kde = ksdensity(data_for_hist, x_kde, 'bandwidth', 0.005); 

hold(ax, 'on'); 
if ~isempty(highlight_ratios)
    for i_cond=1:length(highlight_ratios)
       plot(ax, [highlight_ratios(i_cond), highlight_ratios(i_cond)], ...
                [0, max(kde)*1.1], 'm', 'linew', 4); 
    end
end
fill([x_kde, fliplr(x_kde)], [kde, zeros(size(kde))], col_dark, ...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(x_kde, kde, 'color', col_dark, 'linew', 2)

ax.XLim = [hist_edges(1), hist_edges(end)];ax.XTick = []; 
ax.YAxis.Visible = 'off'; 


% marginal category 2
col_light = [191, 55, 21]/255; 
col_dark = [168, 11, 0]/255; 

ax = pnl(3).select(); 
data_for_hist = reshape(max_ratio_boot(:, par_proto.bound+1:end), [], 1); 
kde = ksdensity(data_for_hist, x_kde, 'bandwidth', 0.005); 

hold(ax, 'on'); 
if ~isempty(highlight_ratios)
    for i_cond=1:length(highlight_ratios)
       plot(ax, [highlight_ratios(i_cond), highlight_ratios(i_cond)], ...
                [0, max(kde)*1.1], 'm', 'linew', 4); 
    end
end
fill([x_kde, fliplr(x_kde)], [kde, zeros(size(kde))], col_dark, ...
     'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(x_kde, kde, 'color', col_dark, 'linew', 2)

ax.XLim = [hist_edges(1), hist_edges(end)];ax.XTick = []; 
ax.YAxis.Visible = 'off'; 


% formating
ax.XTick = ratios_proto; 
ax.XTickLabel = get_ratio_labs(ratios_proto, 'n_decim', 3, 'no_zero', true); 
ax.TickLength = [0,0];  
ax.XTickLabelRotation = -60; 

pnl(1).marginbottom = 4; 
pnl(2).marginbottom = 4; 
pnl(3).marginbottom = 4; 

pnl.margin = [5, 5, 5, 5]; 

if n_prototypes > 20
    pnl.fontsize = 4; 
end

ax = pnl(3).select(); 
for i=1:length(rat_proto_near_zero)
    ax.XTickLabel{idx_proto_with_zero(i)} = ...
            ['\color{red}' ax.XTickLabel{idx_proto_with_zero(i)}];
end


% ---- test against uniform ----

% get a histogram of counts across prototypes 
centres = [1:n_prototypes]; 
centres_ratios = ratios_proto(centres, 1); 
edges = [0.5 : 1 : n_prototypes+0.5]; 
data_for_hist = reshape(max_idx_boot_obs, [], 1); 
counts = histcounts(data_for_hist, edges); 

% select only prototypes between 1:1 and 2:1 (where the stimuli were)
idx_11_21 = dsearchn(ratios_proto(:, 1), stim.ratios([1,end])); 
mask = centres >= idx_11_21(1) & centres <= idx_11_21(2); 

centres_ratios = centres_ratios(mask); 
counts = counts(mask); 

% Chi square goodness of fit against unniform distribution  
N = sum(counts); 
n_bins = length(centres_ratios); 
expected_counts = repmat(N / n_bins, 1, n_bins); 
chi_sq = sum((counts - expected_counts).^2 ./ expected_counts); 
df = n_bins - 1;     
p = chi2cdf(chi_sq, df, 'upper'); 
pnl.title(sprintf('Chisq(%d)=%.2f, p=%.1g', df, chi_sq, p)); 



% ---- on-off stimulus ratio ----

pval = (sum(on_off_index_stim_ratio_boot > 0, 1) + 1) ./ ...
       (size(on_off_index_stim_ratio_boot, 1) + 1); 

pval_adj = min(pval * length(pval), 1); 
   
for i_cond=1:par.n_cond
    
    if pval_adj(i_cond) < 0.05
        
        ax = pnl(2,i_cond).select(); 
        
        ax.Children(end).Color = [167, 75, 184]/255; 
        ax.Children(end).LineWidth = 4; 
        
    end
    
end


% ---- on-off sliding window pre-boundary ----

n_comparisons = size(on_off_index_pre_bound_ratio_boot, 2) + ...
                size(on_off_index_post_bound_ratio_boot, 2);
            
pval = (sum(on_off_index_pre_bound_ratio_boot <= 0, 1) + 1) ./ ...
        (size(on_off_index_pre_bound_ratio_boot, 1) + 1); 

pval_adj = min(pval * n_comparisons, 1); 

prctile(on_off_index_pre_bound_ratio_boot, [2.5, 97.5], 1)'

ax = pnl(1).select(); 
vals = ax.Children(1).YData; 

for i=1:length(pval)
    hold(ax, 'on'); 

    if pval_adj(i) < 0.05
        
        on_lims = par_proto.on_lims_pre;
        off_lims_below = par_proto.off_lims_below_pre;
        off_lims_above = par_proto.off_lims_above_pre;
        
        fill(ax, [on_lims(i,1), on_lims(i,2), on_lims(i,2), on_lims(i,1)], ...
             [0, 0, max(vals)*1.1, max(vals)*1.1], 'k', 'edgecolor', 'none', 'FaceAlpha', 0.2);
        fill(ax, [off_lims_below(i,1), off_lims_below(i,2), off_lims_below(i,2), off_lims_below(i,1)], ...
             [0, 0, max(vals)*1.1, max(vals)*1.1], 'k', 'edgecolor', 'none', 'FaceAlpha', 0.05);
        fill(ax, [off_lims_above(i,1), off_lims_above(i,2), off_lims_above(i,2), off_lims_above(i,1)], ...
             [0, 0, max(vals)*1.1, max(vals)*1.1], 'k', 'edgecolor', 'none', 'FaceAlpha', 0.05);
        
        txt = text(ax, mean([on_lims(i,1), on_lims(i,2)]), max(vals)*1.3, ...
                   sprintf('%.3f\n(%s)', mean([on_lims(i,1), on_lims(i,2)]), ...
                           vec_to_str(on_lims(i,:), 'sep', '-', 'format', '%.2f')), ...
                   'HorizontalAlignment', 'center', 'FontSize', 5); 
         
        plot(ax, [mean([on_lims(i,1), on_lims(i,2)]), mean([on_lims(i,1), on_lims(i,2)])], ...
            [0, max(vals)*1.1], 'k', 'linew', 3); 
         
        ax.YLim = [0, max(vals)*1.4];  
        
    else
        plot(ax, [mean([on_lims(i,1), on_lims(i,2)]), mean([on_lims(i,1), on_lims(i,2)])], ...
            [0, max(vals)*1.1], 'k:', 'linew', 1); 
        
    end
    
end


% ---- on-off sliding window post-boundary ----

n_comparisons = size(on_off_index_pre_bound_ratio_boot, 2) + ...
                size(on_off_index_post_bound_ratio_boot, 2);

pval = (sum(on_off_index_post_bound_ratio_boot <= 0, 1) + 1) ./ ...
        (size(on_off_index_post_bound_ratio_boot, 1) + 1); 

pval_adj = min(pval * n_comparisons, 1); 
    
prctile(on_off_index_post_bound_ratio_boot, [2.5, 97.5], 1)' 

ax = pnl(3).select(); 
vals = ax.Children(1).YData; 

for i=1:length(pval)
    
    hold(ax, 'on'); 

    on_lims = par_proto.on_lims_post;
    off_lims_below = par_proto.off_lims_below_post;
    off_lims_above = par_proto.off_lims_above_post;

    if pval_adj(i) < 0.05
        
        fill(ax, [on_lims(i,1), on_lims(i,2), on_lims(i,2), on_lims(i,1)], ...
             [0, 0, max(vals)*1.1, max(vals)*1.1], 'k', 'edgecolor', 'none', 'FaceAlpha', 0.2);
        fill(ax, [off_lims_below(i,1), off_lims_below(i,2), off_lims_below(i,2), off_lims_below(i,1)], ...
             [0, 0, max(vals)*1.1, max(vals)*1.1], 'k', 'edgecolor', 'none', 'FaceAlpha', 0.05);
        fill(ax, [off_lims_above(i,1), off_lims_above(i,2), off_lims_above(i,2), off_lims_above(i,1)], ...
             [0, 0, max(vals)*1.1, max(vals)*1.1], 'k', 'edgecolor', 'none', 'FaceAlpha', 0.05);
        
        txt = text(ax, mean([on_lims(i,1), on_lims(i,2)]), max(vals)*1.3, ...
                   sprintf('%.3f\n(%s)', mean([on_lims(i,1), on_lims(i,2)]), ...
                           vec_to_str(on_lims(i,:), 'sep', '-', 'format', '%.2f')), ...
                   'HorizontalAlignment', 'center', 'FontSize', 5); 
         
        plot(ax, [mean([on_lims(i,1), on_lims(i,2)]), mean([on_lims(i,1), on_lims(i,2)])], ...
            [0, max(vals)*1.1], 'k', 'linew', 3); 
         
        ax.YLim = [0, max(vals)*1.4];  
        
    else
        plot(ax, [mean([on_lims(i,1), on_lims(i,2)]), mean([on_lims(i,1), on_lims(i,2)])], ...
            [0, max(vals)*1.1], 'k:', 'linew', 1); 
        
    end
    
end


% save distribution and all on-off bootstraps
if do_save
    
    print(fullfile(save_path, [fname_base, '_maxIdxGrandBoot_kde.png']), ...
          '-dpng', '-painters', '-r600', f); 

    print(fullfile(save_path, [fname_base, '_maxIdxGrandBoot_kde.eps']), ...
          '-depsc', '-painters', '-r600', f); 

    save(fullfile(save_path, [fname_base, '_maxIdxGrandBoot_kde.mat']), ...
        'par', 'par_proto', ...
        'feat_proto', 'ratios_proto', 'stim', ...
        'max_idx_boot_obs',  ...
        'on_off_index_stim_ratio_boot', ...
        'on_off_index_pre_bound_ratio_boot', ...
        'on_off_index_post_bound_ratio_boot' ...
        )
    
    close(f); 
end



%%


close all
