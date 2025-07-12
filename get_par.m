function par = get_par(varargin)

parser = inputParser(); 

addParameter(parser, 'experiment', '2ioi'); 

parse(parser, varargin{:}); 

experiment = parser.Results.experiment; 

%% define paths, add libraries

% path to root of the data folder 
data_path = '/Users/tomaslenc/projects_git/Rhototype_2ioi_osf/data'; 

% paths to external libraries 
% ---------------------------

% rnb tools 
rnb_tools_path = '~/projects_git/rnb_tools/src'; 

% letswave 6
lw6_path = '~/projects_git/letswave6'; 

% letswave 
lw7_path = '~/projects_git/letswave7'; 

%% 

% add libraries to matlab path
addpath(genpath(rnb_tools_path)); 
addpath(genpath(lw6_path)); 
addpath(genpath('lib')); 

installBayesFactor

% prepare paths to data subfolders 
raw_path = fullfile(data_path, 'raw'); 
deriv_path = fullfile(data_path, 'derivatives'); 
rsa_path = fullfile(deriv_path, 'rsa'); 
summary_path = fullfile(deriv_path, 'summary'); 

if ~isdir(deriv_path)
    mkdir(deriv_path)
end
if ~isdir(rsa_path)
    mkdir(rsa_path)
end
if ~isdir(summary_path)
    mkdir(summary_path)
end


            
%% parameters

subjects = [6 : 23]; 
% subjects = [15, 20]; 

trial_dur = 22.5; 

cycle_dur = 0.750; 

n_cond = 13; 


% time left before start and after the end of each trial (to allow later
% filtering without edge artifacts)
epoch_buffer_dur = 5; 

% number of channels to use for interpolation 
n_chans_interp = 3; 

filt_cut_low = 0.1; 
filt_cut_high = 64; 
filt_order_low = 4; 
filt_order_high = 4; 

filt_ica_cut_low = 1; 
filt_ica_cut_high = 64; 
filt_ica_order_low = 4; 
filt_ica_order_high = 4; 

decim_factor = 4; 
decim_factor_ica = 4; 

filt_erp_cut_low = 10; 
filt_tap_erp_cut_low = 18; 


%%

roi_front_chans = {
        'F1'
        'Fz'
        'F2'
        'FC1'
        'FCz'
        'FC2'
        'C1'
        'Cz'
        'C2'
    }; 

eeg_ref_chans = {'TP9', 'TP10'}; 

% number of frequency bins for SNR subtraction and zSNR estimation 
snr_bins = [2, 5]; 

spectral_resolution = 1/trial_dur; 
harmonic_resolution = 1/cycle_dur; 

assert((snr_bins(2) * spectral_resolution) < harmonic_resolution, ...
       'The number of SNR bins is too large - there will be overlap to the neighbouring response harmonics'); 

% number of permutations when testing fit of a categorical model at
% individual participant level 
n_perm_rsa_individual = 5000; 
n_perm_rsa_group = 10000; 

% a significant category has to have more pixels (i.e. conditions) than
% this
n_mini_categ_skip = 1; 

%% return structure 

clear varargin

w = whos;
par = []; 
for a = 1:length(w) 
    par.(w(a).name) = eval(w(a).name); 
end










