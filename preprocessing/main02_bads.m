% This script helps to do data cleaning, to visually identify bad channels
% and trials that can be subsequently interpolated or removed respectively.

clear 

par = get_par('experiment', '2ioi'); 

sub = 15; 

%% manual check for bad chans and trials 

wd = pwd; 

% use letswave 7 to open continuous data viewer and check for bad
% trials/channels -> manually write them down 
import_lw(7); 
cd(fullfile(par.deriv_path, 'preproc', sub_num2str(sub))); 
letswave 

warning('set X-range in the viewer to 0 - %.1fs', par.trial_dur); 

%% 

% switch back to letswave 6 when done 
cd(wd); 
import_lw(6);
        
%% generate template table for bad trials and channels 

% This chunk of code will look whether the text file with info about bad
% chanenls and trials exists for the current subject. If not, it will
% generate a new empty one in the proper location and with a proper name.
% Once this is done, open the files manually and enter the bad channels and
% trials you observed during visual inspection in the previous step.

% prepare info about bad channels             
load(fullfile(par.deriv_path, 'bads.mat')); 

idx = find([bads.subject] == sub); 
bad_chans = bads(idx).bad_channels;
bad_epochs = bads(idx).bad_trials;

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% EDIT MANUALLY TO ADD/REMOVE BADS
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% print to console for sanity check 
fprintf('\n\n---------------------------------------------------------\n'); 
bad_chans
bad_epochs
fprintf('---------------------------------------------------------\n'); 

% save back if edited 
bads(idx).bad_channels = bad_chans; 
bads(idx).bad_trials = bad_epochs; 

save(fullfile(par.deriv_path, 'bads.mat'), 'bads'); 

%%

% load EEG data
fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 
fname =  sprintf('sub-%03d_task-EEG_filt_ep', sub); 
[header_ep, data_ep] = CLW_load(fullfile(fpath, fname));

% reject bad trials 
[header_ep, data_ep] = RLW_arrange_epochs(header_ep, data_ep, ...
                            setdiff([1:header_ep.datasize(1)], bad_epochs)); 

% interpolate bad channels 
if ~isempty(bad_chans)
    [header_ep, data_ep] = interpolate_bad_chans(header_ep, data_ep, ...
        bad_chans, par.n_chans_interp); 
    
    bad_chans_str = join(bad_chans, '-');

end

% save
header_ep.name = sprintf('%s_chanint', header_ep.name);
CLW_save(fpath, header_ep, data_ep);


