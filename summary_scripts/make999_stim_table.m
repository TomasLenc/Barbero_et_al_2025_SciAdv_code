% Make table of stimulus inter-onset intraval ratios. 

clear

par = get_par('experiment', '2ioi'); 

save_path = fullfile(par.deriv_path, 'figures', 'tables'); 

if ~isdir(save_path)
    mkdir(save_path)
end

load(fullfile(par.data_path, 'stimuli', 'stim.mat')); 

tbl = table(); 

tbl.condition = [1:13]'; 
tbl.IOI_1 = round(stim.ratios(:, 1) * par.cycle_dur, 3); 
tbl.IOI_2 = round(stim.ratios(:, 2) * par.cycle_dur, 3); 

tbl.IOI_ratio_1 = round(stim.ratios(:, 1), 3); 
tbl.IOI_ratio_2 = round(stim.ratios(:, 2), 3); 

writetable(tbl, fullfile(save_path, 'stim_table.csv'));  