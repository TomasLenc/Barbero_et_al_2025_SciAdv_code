function [r_obs, r_null_dist, p] = ...
    perm_rdm_corr_group_empirical_model(rdms_resp, rdms_model, n_perm, varargin)
% This funciton takes a cell array of response RDMs (corresponding to e.g.
% all subjects in an experiment) and another cell array with model RDMs
% separately for each subject. It performs a group-level permutation test
% by shuffling each response RDM and correlating this to the corresponding
% model to build a null distribution of group-average correlation
% coefficient. 
% 
% Parameters
% ----------
% rdms : cell
%     A cell array containg RDMs. Each element contains the RDM of one
%     participant. 
% model_rdms : cell
%     A cell array containg RDMs. Each element contains the RDM of the
%     model for one participant.
% n_perm : int 
%     How many permutations to run. 
% covariate :  array_like, shape=[n, n], optional
%     RDM that will be regressed out before correlating response RDM and  
%     model model RDM. 
% verbose : bool, optional, default=true
%     If true, info will be printed to the console as the permutations
%     progress. 
% parfor : bool, optional, default=true
%     If true, the permutations will be computed using parallel loop (using
%     all available cores on the machine). If false, a standard for loop
%     will be executed. This is useful if you, for example, want to
%     paralelize the code at a higher level. 
% 
% Returns 
% -------
% r_obs : float 
%     Observed correlation coefficient. 
% r_null_dist : array_like, shape=[1, n_perm]
%     Correlation coefficients of the null distribution. 
% p : float 
%     P-value of the correlation coefficient (test difference from 0)
%     computed as a proportion of more extreme values from the null
%     distribution.
%
parser = inputParser; 

addParameter(parser, 'covariate', []); 
addParameter(parser, 'parfor', true); 
addParameter(parser, 'verbose', true); 

parse(parser, varargin{:}); 

covariate = parser.Results.covariate; 
do_parfor = parser.Results.parfor; 
verbose = parser.Results.verbose; 

%%

assert(length(rdms_resp) == length(rdms_model), ...
    'the number of response RDMs (%d) and model RDMS (%d) doesnt match!', ...
    length(rdms_resp), length(rdms_model))

%%

n_rdms = length(rdms_resp); 

n_cond = size(rdms_resp{1}, 1); 

%% prepare RDM

% get lower-triangular mask 
mask = tril(ones(n_cond, n_cond, 'logical')) & ~diag(ones(1, n_cond, 'logical')) ; 
mask_idx = find(mask); 

% set everything outside of lower triangular to NaN (just in case) 
for i_rdm=1:n_rdms
    rdms_resp{i_rdm}(~mask) = nan; 
    rdms_model{i_rdm}(~mask) = nan; 
end

%% get observed values

r_tmp = nan(1, n_rdms);                      

for i_rdm=1:n_rdms

    % get the correlation for all models 
    r = test_models(rdms_resp{i_rdm}, ...
                    rdms_model{i_rdm}, ...
                    'covariate', covariate, ...
                    'verbose', false);

    % save the correlation 
    r_tmp(i_rdm) = r; 

end

% average across RDMs (subjects) 
r_obs = mean(r_tmp); 

%% do permutations 

r_null_dist = nan(1, n_perm); 

if do_parfor
    parfor i_perm=1:n_perm
        if mod(i_perm, 500) == 0 && verbose
            fprintf('perm %d/%d\n', i_perm, n_perm); 
        end
        r_null_dist(i_perm) = do_one_perm(...
                        mask_idx, rdms_resp, rdms_model, covariate); 
    end
else
    for i_perm=1:n_perm
        if mod(i_perm, 500) == 0 && verbose
            fprintf('perm %d/%d\n', i_perm, n_perm); 
        end
        r_null_dist(i_perm) = do_one_perm(...
                        mask_idx, rdms_resp, rdms_model, covariate); 
    end
end

%% get pvalue from the null distribtuion 

% see: https://stats.stackexchange.com/a/112456
B = sum(r_null_dist > r_obs); 
M = n_perm; 
p = (B+1) ./ (M+1); 


end


function r = do_one_perm(mask_idx, rdms_resp, rdms_model, covariate)

    % randomly shuffle the values on lower triangular of RDM1 
    shuffled_idx = randsample(mask_idx, length(mask_idx), false); 

    r = nan(1, length(rdms_resp));     

    for i_rdm=1:length(rdms_resp)

        shuffled_rdm = nan(size(rdms_resp{i_rdm})); 
        shuffled_rdm(mask_idx) = rdms_resp{i_rdm}(shuffled_idx); 

        % get correlation of each considered model 
        r(i_rdm) = test_models(shuffled_rdm, ...
                               rdms_model{i_rdm}, ...
                               'covariate', covariate, ...
                               'verbose', false);

    end

    % take mean across the group and save it 
    r = mean(r); 

end



