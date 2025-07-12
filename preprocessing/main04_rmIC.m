% This script can be used to visually identify bad IC components and remove
% them from the data

clear 

par = get_par('experiment', '2ioi'); 

sub = 15; 


%%

fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 

%% visual inspecion of ICs 

% open letswave6, and manually inspect the ICs in the
% "sub-$_task-EEG_filt_ep_chanint" dataset, and write down which ones are
% bad

wd = pwd; 
cd(fpath); 
letswave 


%% change dir back to code

cd(wd); 


%% save/load bad ICs 

% prepare info about bad ICs          
load(fullfile(par.deriv_path, 'bads.mat')); 

idx = find([bads.subject] == sub); 
bad_ICs = bads(idx).bad_ICs;

% print to console for sanity check 
fprintf('\n\n---------------------------------------------------------\n'); 
bad_ICs
fprintf('---------------------------------------------------------\n'); 

% save back if edited 
bads(idx).bad_ICs = bad_ICs; 
save(fullfile(par.deriv_path, 'bads.mat'), 'bads'); 

%%

% load the ICA matrix 
clear matrix_ICA
fname = sprintf('sub-%03d_task-EEG_icaMatrix', sub); 
load(fullfile(par.deriv_path, fname));

% load data 
fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint', sub); 
[header_ep, data_ep] = CLW_load(fullfile(fpath, fname));

% remove bad ICs from data 
n_ic = size(matrix_ICA.ica_mm, 2); 
ic2keep = setdiff([1 : n_ic], bad_ICs); 

[header_ep, data_ep] = RLW_ICA_apply_filter(header_ep, data_ep, ...
                                      matrix_ICA.ica_mm, ...
                                      matrix_ICA.ica_um, ...
                                      ic2keep);     

% save clean data 
header_ep.name = [header_ep.name, sprintf('_icfilt(%s)', ...
                                          vec_to_str(bad_ICs, 'sep',','))];

CLW_save(fpath, header_ep, data_ep); 


