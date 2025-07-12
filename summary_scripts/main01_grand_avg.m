% Compute grand averages for different kinds of EEG and tapping responses. 

clear

par = get_par('experiment', '2ioi'); 

subjects = par.subjects;  

do_eeg_cont     = true; 
do_eeg_erp      = true; 
do_tap_force    = true; 
do_tap_impulse  = true; 
do_tap_erp      = true; 
do_tap_ratios   = true; 

%%

data_eeg_all = cell(1, length(subjects)); 
data_tap_all = cell(1, length(subjects)); 
data_tap_imp_all = cell(1, length(subjects)); 
data_tap_ratios_all = cell(1, length(subjects)); 

header_erp_all = cell(1, length(subjects)); 
data_erp_all = cell(1, length(subjects)); 

header_tap_env_erp_all = cell(1, length(subjects)); 
data_tap_env_erp_all = cell(1, length(subjects)); 

for i_sub=1:length(subjects)

    sub = subjects(i_sub); 

    sub_str = sub_num2str(sub)

    if do_eeg_cont
        % low-frequency response 
        % ----------------------

        % load data 
        fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
        fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt*.lw6', sub); 
        d = dir(fullfile(fpath, fname)); 

        [header_eeg_ep, data_eeg_ep] = CLW_load(fullfile(d.folder, d.name));

        % rereference 
        [header_eeg_ep, data_eeg_ep] = RLW_rereference(header_eeg_ep, data_eeg_ep, ...
                                        'apply_list', {header_eeg_ep.chanlocs.labels}, ...
                                        'reference_list', par.eeg_ref_chans); 

        % don't concatenate epochs
        [header_eeg_avg, data_eeg_avg] = average_trials_per_cond(header_eeg_ep, data_eeg_ep, ...
                                    'concatenate_epochs', false); 

        data_eeg_all{i_sub} = data_eeg_avg; 

    end
        
    if do_eeg_erp
        
        % ERP response 
        % ----------------------

        % load data 
        fpath = fullfile(par.deriv_path, 'erp', sub_str); 
        fname =  sprintf('sub-%03d_response-erp_lpf-%.0fHz', ...
                         sub, ...
                         par.filt_erp_cut_low); 

        [header_erp, data_erp] = CLW_load(fullfile(fpath, fname));

        [header_erp.events.sub] = deal(sub);                 

        header_erp_all{i_sub} = header_erp; 
        data_erp_all{i_sub} = data_erp; 
        
    end
    

    if do_tap_force
    
        % tapping response 
        % -----------------

        tap_response_type = 'force'; 
        
        % load data 
        fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
        fname =  sprintf('sub-%03d_task-TAP_response-%s.lw6', ...
                         sub, tap_response_type); 
        d = dir(fullfile(fpath, fname)); 

        [header_tap_env, data_tap_env] = CLW_load(fullfile(d.folder, d.name));

        % segment all trials
        codes = unique({header_tap_env.events.code}); 
        [header_tap_env_ep, data_tap_env_ep] = segment_safe(...
                                                header_tap_env, ...
                                                data_tap_env, ...
                                                codes,...
                                                'x_start', 0,  ...
                                                'x_duration', par.trial_dur); 

        % don't concatenate epochs
        [header_tap_avg, data_tap_avg] = average_trials_per_cond(...
                                    header_tap_env_ep, data_tap_env_ep, ...
                                    'concatenate_epochs', false); 

        data_tap_all{i_sub} = data_tap_avg; 
        
    end
    
    
    
   if do_tap_impulse
    
        % tapping impulse timeseries 
        % -----------------
        
        % load data 
        fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
        fname =  sprintf('sub-%03d_task-TAP_response-impulse.lw6', ...
                         sub); 
        d = dir(fullfile(fpath, fname)); 

        [header_tap_imp, data_tap_imp] = CLW_load(fullfile(d.folder, d.name));

        % segment all trials
        codes = unique({header_tap_imp.events.code}); 
        [header_tap_imp_ep, data_tap_imp_ep] = segment_safe(...
                                                header_tap_imp, ...
                                                data_tap_imp, ...
                                                codes,...
                                                'x_start', 0,  ...
                                                'x_duration', par.trial_dur); 

        % don't concatenate epochs
        [header_tap_imp_avg, data_tap_imp_avg] = average_trials_per_cond(...
                                    header_tap_imp_ep, data_tap_imp_ep, ...
                                    'concatenate_epochs', false); 

        data_tap_imp_all{i_sub} = data_tap_imp_avg; 
        
    end    
    
    
    if do_tap_erp
        
        % tap envelope ERP response 
        % ----------------------

        % load data 
        fpath = fullfile(par.deriv_path, 'erp', sub_str); 
        fname =  sprintf('sub-%03d_response-tap-force-erp_lpf-%.0fHz', ...
                         sub, ...
                         par.filt_tap_erp_cut_low); 

        [header_tap_env_erp, data_tap_env_erp] = CLW_load(fullfile(fpath, fname));

        [header_tap_env_erp.events.sub] = deal(sub);                 

        header_tap_env_erp_all{i_sub} = header_tap_env_erp; 
        data_tap_env_erp_all{i_sub} = data_tap_env_erp; 
        
    end
    
    
    if do_tap_ratios
        
        % tap onset ratios
        % ----------------------
        
        fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
        fname =  sprintf('sub-%03d_response-tap-ratios.mat', sub); 
        data_tap_ratios = load(fullfile(fpath, fname)); 
        
        data_tap_ratios_all{i_sub} = data_tap_ratios.tap_ratios; 
        
    end
    
end


%% EEG ERP 

if do_eeg_erp

    data_erp = cat(1, data_erp_all{:}); 

    events = cellfun(@(x) x.events, header_erp_all, 'uni', 0); 
    header_erp = header_erp_all{1}; 
    header_erp.events = cat(2, events{:}); 
    header_erp.datasize = size(data_erp); 

    [header_erp_mean, data_erp_mean] = average_trials_per_cond(...
                                                header_erp, data_erp); 

    save_path = fullfile(par.deriv_path, 'erp', 'sub-grand'); 
    mkdir(save_path); 

    fname = sprintf('sub-grand_response-erp_lpf-%.0fHz', par.filt_erp_cut_low);     
    header_erp_mean.name = fname; 
    CLW_save(save_path, header_erp_mean, data_erp_mean); 

    fname = sprintf('sub-grand_response-erp_lpf-%.0fHz_ind', par.filt_erp_cut_low);     
    header_erp.name = fname; 
    CLW_save(save_path, header_erp, data_erp); 

    % %sanity check (pretty)
    % load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'))
    % plot_erp_overlay(squeeze(mean(data_erp_mean, 2)), 1/header_erp.xstep, ...
    %                  'ratios', stim.ratios)
    
    
end

%% EEG FFT

if do_eeg_cont

    header_lf = header_eeg_avg; 
    data_lf = mean(cat(7, data_eeg_all{:}), 7); 

    % re-segment all trials
    codes = unique({header_lf.events.code}); 

    [header_lf, data_lf] = segment_safe(header_lf, data_lf, ...
                                              codes,...
                                              'x_start', 0, ...
                                              'x_duration', par.trial_dur); 

    % get FFTs and save them                     
    save_path = fullfile(par.deriv_path, 'fft', 'sub-grand'); 
    mkdir(save_path); 

    fname = sprintf('sub-grand_response-eeg_concatEp-false'); 

    get_and_save_ffts(header_lf, data_lf, par.snr_bins, ...
                      'fname', fname, 'save_path', save_path); 



end


%% TAP FFT

if do_tap_force
    
    header_tap = header_tap_avg; 

    data_tap = mean(cat(7, data_tap_all{:}), 7); 

    save_path = fullfile(par.deriv_path, 'fft', 'sub-grand'); 
    mkdir(save_path); 

    fname = sprintf('sub-grand_response-tap-force_concatEp-false'); 

    get_and_save_ffts(header_tap, data_tap, par.snr_bins, ...
                      'fname', fname, 'save_path', save_path, 'fmax', 100); 

end


%% TAP IMPULSE FFT

if do_tap_impulse
    
    header_tap = header_tap_imp_avg; 

    data_tap = mean(cat(7, data_tap_imp_all{:}), 7); 

    save_path = fullfile(par.deriv_path, 'fft', 'sub-grand'); 
    mkdir(save_path); 

    fname = sprintf('sub-grand_response-tap-impulse_concatEp-false'); 

    get_and_save_ffts(header_tap, data_tap, par.snr_bins, ...
                      'fname', fname, 'save_path', save_path); 

end

%% TAP FORCE ERP 

if do_tap_erp

    data_tap_env_erp = cat(1, data_tap_env_erp_all{:}); 

    events = cellfun(@(x) x.events, header_tap_env_erp_all, 'uni', 0); 
    
    header_tap_env_erp = header_tap_env_erp_all{1}; 
    header_tap_env_erp.events = cat(2, events{:}); 
    header_tap_env_erp.datasize = size(data_tap_env_erp); 

    [header_tap_env_erp_mean, data_tap_env_erp_mean] = average_trials_per_cond(...
                                    header_tap_env_erp, data_tap_env_erp); 

    save_path = fullfile(par.deriv_path, 'erp', 'sub-grand'); 
    mkdir(save_path); 

    fname = sprintf('sub-grand_response-tap-force-erp_lpf-%.0fHz', par.filt_tap_erp_cut_low);     
    header_tap_env_erp_mean.name = fname; 
    CLW_save(save_path, header_tap_env_erp_mean, data_tap_env_erp_mean); 

    fname = sprintf('sub-grand_response-tap-force-erp_lpf-%.0fHz_ind', par.filt_tap_erp_cut_low);     
    header_tap_env_erp.name = fname; 
    CLW_save(save_path, header_tap_env_erp, data_tap_env_erp); 

 
end


%% TAP RATIOS

if do_tap_ratios
    
    tap_ratios = mean(cat(3, data_tap_ratios_all{:}), 3); 
    
    save_path = fullfile(par.deriv_path, 'preproc', 'sub-grand'); 
    mkdir(save_path); 

    fname =  sprintf('sub-grand_response-tap-ratios.mat'); 
    
    save(fullfile(save_path, fname), 'tap_ratios'); 
    
end




