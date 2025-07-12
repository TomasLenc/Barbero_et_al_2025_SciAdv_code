% Plot RSM using stimulus inter-onset interval ratios

clear 

par = get_par('experiment', '2ioi'); 

flip_rdm = false; 

%% load

fpath = fullfile(par.data_path, 'stimuli'); 
load(fullfile(fpath, 'stim.mat')); 

%% interval ratios

rdm = get_rdm_from_feat(stim.ratios(:,1)', 'method', 'L1');

f = plot_rdm(rdm, 'ratios', stim.ratios, 'flip', flip_rdm); 

save_path_onsets = fullfile(par.deriv_path, 'results', 'response-stim-ratios'); 
mkdir(save_path_onsets); 

fname = sprintf('response-stim_feat-ratios_method-pearson_rdm'); 

save(fullfile(save_path_onsets, [fname, '.mat']), 'rdm'); 
print(f, '-dpng', '-painters', '-r300', fullfile(save_path_onsets, fname))

close(f)
