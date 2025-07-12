% Save tapping force as a separate file. Use tap onset times to prepare
% continuous onset vectors.

clear

par = get_par('experiment', '2ioi'); 

sub = 15; 

%% copy the force siganl from raw

sub_str = sub_num2str(sub); 

fpath = fullfile(par.raw_path); 

fname = sprintf('%s_task-TAP.lw6', sub_str); 

[header, data] = CLW_load(fullfile(fpath, fname)); 

data = data(:, 1, 1, 1, 1, :); 
header.datasize = size(data); 
header.chanlocs = header.chanlocs(1); 

header.name = sprintf('sub-%03d_task-TAP_response-force', sub);

CLW_save(fullfile(par.deriv_path, 'preproc', sub_str), header, data); 


%% crate continous tap onset signals 

% load data
fpath = fullfile(par.deriv_path, 'preproc', sub_str); 

fname = sprintf('sub-%03d_response-tap-ratios.mat', sub); 

data_onsets = load(fullfile(fpath, fname))

fs = 256; 

N = round(fs * par.trial_dur); 


c_event = 1; 
events = []; 
data_all = {par.n_cond, 1}; 

for i_cond=1:par.n_cond
    
    data_cond = zeros(length(data_onsets.tap_onset_times_all{i_cond}), 1, 1, 1, 1, N); 

    for i_ep=1:length(data_onsets.tap_onset_times_all{i_cond})
        
        tap_onset_times = data_onsets.tap_onset_times_all{i_cond}{i_ep}; 
        
        idx = round(tap_onset_times * fs) + 1; 
        
        % there may be a case where the very last sample contained a tap,
        % lets remove that tap 
        idx(idx == N+1) = []; 
                    
        if any(idx > N)
            error('some tap onsets beyond trial duration???')
        end
        
        data_cond(i_ep, 1, 1, 1, 1, idx) = 1; 
        
        events(c_event).code = num2str(i_cond); 
        events(c_event).latency = 0; 
        events(c_event).epoch = c_event; 

        c_event = c_event+1; 
    end
    
    data_all{i_cond} = data_cond; 

end

data = cat(1, data_all{:}); 
  
header = make_lw_header(...
                    'data', data, ...
                    'fs', fs, ...
                    'events', events); 
                    
header.chanlocs.labels = 'onsets'; 

% save
fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
mkdir(fpath); 

header.name = sprintf('sub-%03d_task-TAP_response-impulse', sub);
CLW_save(fpath, header, data); 

