% Generate figures of stimulus FFT and the corresponding RSMs. 

clear  

par = get_par('experiment', '2ioi'); % 

% flip the order of conditions for plotting
flip_rdm = false; 

% select harmonics for the RDM
harmonics = [1 : 6]; 

%% prepare paths

harm_str = vec_to_str(harmonics, 'format','%d', 'sep','-'); 

frex = 1 / par.cycle_dur * harmonics; 

save_path = fullfile(par.deriv_path, 'results', ...
                     'response-stim', ...
                     sprintf('harm-%s', harm_str)); 
                 
mkdir(save_path); 


%% load

fpath = fullfile(par.data_path, 'stimuli'); 
load(fullfile(fpath, 'stim.mat')); 


%% plot time-domain and magnitude spectrum

n_cond = length(stim.s); 

f = figure('color', 'white', 'position', [631 71 600 930]); 
pnl = panel(f); 
pnl.pack('h', 2); 
pnl(1).pack('v', n_cond); 
pnl(2).pack('v', n_cond); 

pnl.de.margin = 0;

y_min = Inf; 
y_max = -Inf;

N = length(stim.s{1}); 
fs = stim.fs; 
freq = [0 : N-1] / N * fs;   

t = [0 : N-1] / fs; 

fmin = 0.1; 
fmax = 15; 

for i_cond=1:n_cond
               
    frex_idx = round(frex / fs * N) + 1; 
    min_freq_idx = round(fmin / fs * N) + 1; 
    max_freq_idx = round(fmax / fs * N) + 1; 
        
    s = stim.s{i_cond}; 
    env = abs(hilbert(s)); 
    
    mX_to_plot = abs(fft(env)) ./ N * 2; 
    mX_to_plot(1) = 0; 
    
    ax = pnl(1, i_cond).select(); 
    
    plot(t, s, 'linew', 0.6); 
    hold(ax, 'on'); 
    plot(t, env, 'r', 'linew', 0.5); 
    
    ax.XLim = [0, stim.cycle_dur * 2]; 
    ax.YTick = [];     
    
    if i_cond ~= n_cond
        ax.XTick = []; 
        ax.YTick = []; 
    end
    
    ax = pnl(2, i_cond).select(); 

    plot(freq(min_freq_idx:max_freq_idx), ...
         mX_to_plot(min_freq_idx:max_freq_idx), ...
         'linew', 1); 
     
    hold on 
        
    plot(freq(frex_idx), mX_to_plot(frex_idx), 'ro'); 
    
    if i_cond ~= n_cond
        ax.XTick = []; 
        ax.YTick = []; 
        ax.XAxis.Visible = 'off'; 
    end
        
    y_min = min(y_min, min(mX_to_plot(min_freq_idx:max_freq_idx))); 
    y_max = max(y_max, max(mX_to_plot(min_freq_idx:max_freq_idx))); 
    
end

for i_cond=1:n_cond
    ax = pnl(2, i_cond).select(); 
    ax.YLim = [y_min, y_max]; 
    ax.Title.Position(2) = ax.Title.Position(2) - ax.Title.Position(2) * 0.3; 
end

ax.YTick = [y_min, y_max]; 
ax.YTickLabel = [round(y_min * 100)/100, round(y_max * 100)/100]; 
ax.TickDir = 'out'; 

pnl.de.margintop = 8;

pnl.de.marginleft = 20; 

pnl.marginbottom = 10; 

pnl.margintop = 10; 

pnl.fontsize = 12; 


fname = sprintf('response-stim-env_time_mX'); 

print(f, '-dpng', '-painters', '-r300', fullfile(save_path, fname))

close(f); 

%% RSA

% complex spectrum representation
stim_env = cellfun(@(x) abs(hilbert(x)), stim.s, 'uni', 0); 
[X_stim, freq_stim] = get_X(stim_env, stim.fs); 


feat_stim = get_feat_from_X(X_stim, freq_stim, frex); 
rdm = get_rdm_from_feat(feat_stim, 'method', 'pearson'); 
f = plot_rdm(rdm, 'ratios', stim.ratios, 'flip', flip_rdm); 

fname = sprintf('response-stim-env_feat-X_method-pearson_rdm'); 

save(fullfile(save_path, [fname, '.mat']), 'rdm'); 

print(f, '-dpng', '-painters', '-r300', fullfile(save_path, fname))

close(f)

