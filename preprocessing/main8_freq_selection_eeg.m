% Analysis of EEG responses at individual frequencies of interest to
% determine the frequency range of the response. 

clear 

%%

experiment = '2ioi'; 

response = 'eeg'; 

roi_name = 'all';         
avg_chans = false; 

% roi_name = 'front';         
% avg_chans = true; 

concat_ep = false; 

x_lims = [0.5, 30]; 

max_harm = floor(x_lims(2) / (1/0.750)); 


%%

par = get_par('experiment', experiment); 

subjects = num2cell(par.subjects); 

snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2)); 


%%


subdirs = fullfile(sprintf('response-%s', response), ...
                         sprintf('harm-all'), ...
                         sprintf('concatEp-%s_snr-%s', jsonencode(concat_ep), snr_str), ...
                         sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)) ...
                         ); 

load_path = fullfile(par.rsa_path, subdirs); 

save_path = fullfile(par.summary_path, subdirs); 
mkdir(save_path); 

fname_save_base = fullfile(save_path, ...
        sprintf('response-%s_feat-mX_harm-all_concatEp-%s_snr-%s_roi-%s_avgChans-%s', ...
                         response, ...
                         jsonencode(concat_ep), ...
                         snr_str, ...
                         roi_name, ...
                         jsonencode(avg_chans) ...
                         )); 

%% prepare grand average mX  

data_all = cell(1, length(subjects)); 

for i_sub=1:length(subjects)
    sub = subjects{i_sub}; 
    sub_str = sub_num2str(sub); 
    % load data 
    load_path = fullfile(par.deriv_path, 'fft', sub_str); 
    fname = sprintf('sub-%03d_response-eeg_concatEp-false_snr-none_mX', sub); 
    [header, data] = CLW_load(fullfile(load_path, fname));
    % remove mastoids 
    mask_mast = ismember({header.chanlocs.labels}, {'TP9', 'TP10'}); 
    [header, data] = RLW_arrange_channels(header, data, ...
        {header.chanlocs(find(~mask_mast)).labels}); 
    % append the data
    data_all{i_sub} = data; 
end

chanlocs = header.chanlocs(~mask_mast); 

% average across all subjects and conditions 
mX_all_chan = squeeze(mean(cat(1, data_all{:}), 1)); 

% average across all channels 
mX = mean(mX_all_chan, 1); 

% subtract SNR 
freq = [0 : header.datasize(end)-1] * header.xstep; 



%% z-snr by 1 frequency

harm = [1 : max_harm]; 
frex = harm * 1/par.cycle_dur; 
z_snr = nan(1, length(frex)); 
for i_f=1:length(frex)
    z_snr(i_f) = get_z_snr(mX, freq, frex(i_f), ...
                               par.snr_bins(1), par.snr_bins(2)); 
end
[f, ax] = plot_val_over_f(frex, z_snr); 
hold on 
z_thr = norminv(1 - 1e-3); 
plot([0, 60], [z_thr, z_thr], '--', 'Color', [0.5, 0.5, 0.5], 'linew', 2); 
mask_signif = z_snr > z_thr; 
plot(frex(mask_signif), z_snr(mask_signif) + 0.4*mean(z_snr), '*')
txt = text(x_lims(2), max(z_snr)*0.8, sprintf('zSNR thr = %.2f', z_thr)); 
txt.HorizontalAlignment = 'right'; 
txt.FontSize = 12; 
ax.XLim = x_lims; 

ax.YLim = [0, 100]; 

saveas(f, [fname_save_base, '_zSnrIndBy1Frex.svg']); 

tbl = cell2table(num2cell([harm', frex', z_snr', mask_signif']),...
                 'VariableNames', {'harmonic', 'frex', 'z_snr', 'signif'}); 
writetable(tbl, [fname_save_base, '_zSnrIndBy1Frex.csv']); 



%% topo colormap 

addpath(genpath('lib/cmaps')); 
cmap = viridis(256); 

f = figure('color', 'white', 'pos', [476 97 276 609]); 
ax = axes(f); 
for i=1:size(cmap, 1)
    fill([0, 1, 1, 0], [0, 0, 1, 1]+i*0.8, cmap(i, :), 'LineStyle', 'none'); 
    hold on
end
ax.Visible = 'off'; 

print(f, '-dpng', '-painters', '-r600', [fname_save_base, '_topoCbar.png'])
close(f); 

%% plot empty topo head with frontal ROI selection 

f = figure('pos',[620 882 187 165], 'color', 'white'); 

topoplot(get_chan_idx(chanlocs, par.roi_front_chans), chanlocs, ...
         'style', 'blank', ...
         'whitebk', 'on', ...
         'electrodes', 'off');


saveas(f, [fname_save_base, '_topoRoiSelFront.svg']); 

%% plot empty topo head with all channels selection

f = figure('pos',[620 882 187 165], 'color', 'white'); 

topoplot(ones(1,64), chanlocs, ...
         'style', 'blank', ...
         'whitebk', 'on', ...
         'electrodes', 'off');


saveas(f, [fname_save_base, '_topoRoiSelAll.svg']); 




%%

close all


