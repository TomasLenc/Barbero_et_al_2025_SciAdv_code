function [r_obs, r_null_dist, p, models, models_info] = ...
    perm_rdm_corr_group_fix_model(rdms, n_categ, n_perm, varargin)
% This funciton takes a cell array of RDMs (corresponding to e.g. all
% subjects in an experiment) and performs a group-level permutation test.
% It finds the correlation of each possible categorical model with each RDM
% separately, and takes the average correlation coefficient across
% the group. Then, the function generates a null distribution of
% group-average correlation coefficients, by permuting the RDMs of all
% participants (within-participant permutation, but the same random
% shuffling is applied to all input RDMs at one partiuclar permutation
% iteration).
% 
% Parameters
% ----------
% rdms : cell
%     A cell array containg RDMs. Each element contains the RDM of one
%     participant. 
% n_categ : int {2, 3}
%     Number of categories to fit in the model. 
% n_perm : int 
%     How many permutations to run. 
% n_mini_categ_skip : int, optional, default=0
%     Categories with <= n_mini_categ_skip pixels will be excluded. 
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

addParameter(parser, 'n_mini_categ_skip', 0); 
addParameter(parser, 'p_adj_method', 'bonf'); 
addParameter(parser, 'covariate', []); 
addParameter(parser, 'parfor', true); 
addParameter(parser, 'verbose', true); 

parse(parser, varargin{:}); 

n_mini_categ_skip = parser.Results.n_mini_categ_skip; 
p_adj_method = parser.Results.p_adj_method; 
covariate = parser.Results.covariate; 
do_parfor = parser.Results.parfor; 
verbose = parser.Results.verbose; 

%%

n_rdms = length(rdms); 

n_cond = size(rdms{1}, 1); 

%% prepare RDM

% get lower-triangular mask 
mask = tril(ones(n_cond, n_cond, 'logical')) & ~diag(ones(1, n_cond, 'logical')) ; 
mask_idx = find(mask); 

% set everything outside of lower triangular to NaN (just in case) 
for i_rdm=1:n_rdms
    rdms{i_rdm}(~mask) = nan; 
end

%% generate all possible categorical models 

[models, models_info] = generate_categ_models(n_cond, n_categ, ...
                                'n_mini_categ_skip', n_mini_categ_skip); 

                            
r_obs = nan(1, length(models)); 
p = nan(1, length(models)); 
r_null_dist = nan(n_perm, length(models)); 

%% get observed values

for i_model=1:length(models)
    
    r_tmp = nan(1, n_rdms);                      

    for i_rdm=1:n_rdms

        % get the correlation for all models 
        r = test_models(rdms{i_rdm}, ...
                        models(i_model), ...
                        'covariate', covariate, ...
                        'verbose', false);

        % save the correlation of the best fitting model 
        r_tmp(i_rdm) = max(r); 

    end

    % average across RDMs (subjects) 
    r_obs(i_model) = mean(r_tmp); 
    
end

%% do permutations 

if do_parfor
    parfor i_perm=1:n_perm
        if mod(i_perm, 500) == 0 && verbose
            fprintf('perm %d/%d\n', i_perm, n_perm); 
        end
        r_null_dist(i_perm, :) = do_one_perm(...
                        mask_idx, rdms, models, covariate); 
    end
else
    for i_perm=1:n_perm
        if mod(i_perm, 500) == 0 && verbose
            fprintf('perm %d/%d\n', i_perm, n_perm); 
        end
        r_null_dist(i_perm, :) = do_one_perm(...
                        mask_idx, rdms, models, covariate); 
    end
end

%% get pvalue from the null distribtuion 

% see: https://stats.stackexchange.com/a/112456
B = sum(r_null_dist > r_obs, 1); 
M = n_perm; 
p = (B+1) ./ (M+1); 

if strcmpi(p_adj_method, 'bonf')
    p = min(p .* length(models), 1); 
end



end


function r = do_one_perm(mask_idx, rdms, models, covariate)

    % randomly shuffle the values on lower triangular of RDM1 
    shuffled_idx = randsample(mask_idx, length(mask_idx), false); 

    r = nan(length(rdms), length(models));     

    for i_rdm=1:length(rdms)
        for i_model=1:length(models)

            shuffled_rdm = nan(size(rdms{i_rdm})); 
            shuffled_rdm(mask_idx) = rdms{i_rdm}(shuffled_idx); 

            % get correlation of each considered model 
            r(i_rdm, i_model) = test_models(shuffled_rdm, ...
                           models(i_model), ...
                           'covariate', covariate, ...
                           'verbose', false);

        end
    end

    % take mean across the group and save it 
    r = mean(r, 1); 

end



