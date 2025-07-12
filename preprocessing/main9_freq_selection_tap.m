% Analysis of tapping responses at individual frequencies of interest to
% determine the frequency range of the response. 

clear 

%%

experiment = '2ioi'; 

response = 'tap-force'; % tap-force, tap-impulse

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
                         sprintf('concatEp-%s_snr-%s', jsonencode(concat_ep), snr_str) ...
                         ); 

load_path = fullfile(par.rsa_path, subdirs); 

save_path = fullfile(par.summary_path, subdirs); 
mkdir(save_path); 

fname_save_base = fullfile(save_path, ...
        sprintf('response-%s_feat-mX_harm-all_concatEp-%s_snr-%s', ...
                         response, ...
                         jsonencode(concat_ep), ...
                         snr_str)); 




%% prepare grand average mX  

data_all = cell(1, length(subjects)); 

for i_sub=1:length(subjects)
    sub = subjects{i_sub}; 
    sub_str = sub_num2str(sub); 
    % load data 
    fpath = fullfile(par.deriv_path, 'fft', sub_str); 
    fname = sprintf('sub-%03d_response-%s_concatEp-false_snr-none_mX', sub, response); 
    [header, data] = CLW_load(fullfile(fpath, fname));
    % append the data
    data_all{i_sub} = data; 
end

% average across all channels 
mX = squeeze(mean(cat(1, data_all{:}), 1))'; 

% subtract SNR 
freq = [0 : header.datasize(end)-1] * header.xstep; 

mX_snr = subtract_noise_bins(mX, par.snr_bins(1), par.snr_bins(2)); 


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
xlabel('freq (Hz)')
ylabel('zSNR')

ax.YLim = [0, 100]; 

saveas(f, [fname_save_base, '_zSnrIndBy1Frex.svg']); 

tbl = cell2table(num2cell([harm', frex', z_snr', mask_signif']),...
                 'VariableNames', {'harmonic', 'frex', 'z_snr', 'signif'}); 
writetable(tbl, [fname_save_base, '_zSnrIndBy1Frex.csv']); 

%% 

close all
