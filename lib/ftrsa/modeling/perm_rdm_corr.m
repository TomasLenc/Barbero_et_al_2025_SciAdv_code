function [r_obs, p_perm, r_null_dist] = perm_rdm_corr(rdm_resp, rdm_model, n_perm, varargin)
% This function will take two RDMs and do a permutation test, testing
% whether the two matrices are different or not.
%
% Parameters
% ----------
% rdm : array_like, shape=[n, n]
%     Representational dissimiliary matrix with n dimensions. 
% n_perm : int
%     How many permutations to run. 
% 
% Returns 
% -------
% r_obs : float 
%     Observed correlation coefficient between the RDM and best-fitting
%     categorical model. 
% p_perm : float 
%     P-value of the correlation coefficient computed as a proportion of
%     more extreme values in the null distribution.
% r_null_dist : array_like, shape=[1, n_perm]
%     Correlation coefficients of the null distribution. 

parser = inputParser; 

addParameter(parser, 'covariate', []); 
addParameter(parser, 'parfor', true); 
addParameter(parser, 'verbose', true); 

parse(parser, varargin{:}); 

covariate = parser.Results.covariate; 
do_parfor = parser.Results.parfor; 
verbose = parser.Results.verbose; 


%% input checking

% Sanity check - there should be no nans. 
assert(all(~isnan(rdm_resp(:))), ...
    'find_test_best_categ_model:nansInRDM', ...
    'There are NaNs in the RDM!'); 

% It can happen that the data and model matrix could be incompatible, in
% the sense that one is a similarity matrix (e.g. estimated using
% correlation across conditions), while the other is a dissimmilarity
% matrix (e.g. estimated using L2 norm across conditions). Since this
% funciton tests models that are based on similarity, the RDM of the data
% should have ones on the diagonal. 
assert(all(diag(rdm_resp) == 1), ...
    'find_test_best_categ_model:dataNotSimilarityMatrix', ...
    'The input seems to be a dissimilarity matrix, but you need to provide a similarity matrix!'); 

%% prepare RDM

n_cond = size(rdm_resp, 1); 

% get lower-triangular mask 
mask = tril(ones(n_cond, n_cond, 'logical')) & ~diag(ones(1, n_cond, 'logical')) ; 
mask_idx = find(mask); 

% set everything outside of lower triangular to NaN (just in case) 
rdm_resp(~mask) = nan; 
rdm_model(~mask) = nan; 

%% get observed values
                                                    
r_obs = get_rdm_model_corr(rdm_resp, rdm_model, ...
                           'covariate', covariate, ...
                           'verbose', false);

%% do permutations to test significance of the best model

r_null_dist = nan(1, n_perm); 

if do_parfor
    parfor i_perm=1:n_perm
        if mod(i_perm, 500) == 0 && verbose
            fprintf('perm %d/%d\n', i_perm, n_perm); 
        end
        r_null_dist(i_perm) = do_one_perm(...
                        mask_idx, rdm_resp, rdm_model, covariate); 
    end
else
    for i_perm=1:n_perm
        if mod(i_perm, 500) == 0 && verbose
            fprintf('perm %d/%d\n', i_perm, n_perm); 
        end
        r_null_dist(i_perm) = do_one_perm(...
                        mask_idx, rdm_resp, rdm_model, covariate); 
    end
end

%% get pvalue from the null distribtuion 

% see: https://stats.stackexchange.com/a/112456
B = sum(r_null_dist > r_obs); 
M = n_perm; 
p_perm = (B+1) / (M+1); 


end


function r = do_one_perm(mask_idx, rdm_resp, rdm_model, covariate)

    % randomly shuffle the values on lower triangular of RDM1 
    shuffled_idx = randsample(mask_idx, length(mask_idx), false); 

    shuffled_rdm = nan(size(rdm_resp)); 
    shuffled_rdm(mask_idx) = rdm_resp(shuffled_idx); 

    % get correlation of each considered model 
    r = test_models(shuffled_rdm, ...
                    rdm_model, ...
                    'covariate', covariate, ...
                    'verbose', false);

end




