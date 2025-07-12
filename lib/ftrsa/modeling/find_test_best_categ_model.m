function [r_obs, p_perm, r_null_dist, best_model, best_model_info, all_models_info] = ...
                        find_test_best_categ_model(rdm, n_categ, n_perm, varargin)
% Top-level function that takes an RDM and finds the best-fitting
% categorical model from a range of categorical models that differ in the
% position of the category boundary. Once the best model is found, its
% significance is estimated using permutation testing whereby the category
% labels in the RDM are randomly shuffled, and the permuted RDM is
% correlated with the best-fitting categorical model to build a null
% distribution of correlation coefficients.  
% 
% Parameters
% ----------
% rdm : array_like, shape=[n, n]
%     Representational dissimiliary matrix with n dimensions. 
% n_categ : array of int {2, 3, [2, 3]}
%     Number of categories to fit in the model. It can be all 2-category
%     models, all 3-category models, or all 2 AND 3 category models. 
% n_perm : int
%     How many permutations to run. 
% p_adj_method : str, {'none', 'bonf'}, optional, default='bonf'
%     Correction for multiple comparisons. Default is Bonferroni
%     correction, i.e. multiplying the observed permutation pvalue by the
%     number of tested models.
% n_mini_categ_skip : int, optional, default=0
%     Models that include categories with <= `n_mini_categ_skip` pixels
%     will be excluded.
% covariate :  array_like, shape=[n, n], optional
%     DSM that will be regressed out before correlating rdm and categorical 
%     model. 
% verbose : bool, optional, default=true
%     If true, info will be printed to the console as the permutations
%     progress. 
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
% best_model : array_like, shape=[n, n]
%     RDM of the best-fitting categorical model. 
% best_model_info : table, shape=[1, k]
%     One row of table (k columns) which describes the parameters of the
%     winnig model, and its fit. 
% all_models_info : table, shape=[n_models, k]
%     Table (k columns) which describes the parameters of each tested
%     categorical model. Each row is a model, and the rows are ordered from
%     best-fitting to worst-fitting model based on the observed correlation
%     coefficient. 
% 

parser = inputParser; 

addParameter(parser, 'n_mini_categ_skip', 0); 
addParameter(parser, 'p_adj_method', 'bonf'); 
addParameter(parser, 'covariate', []); 
addParameter(parser, 'verbose', true); 

parse(parser, varargin{:}); 

n_mini_categ_skip = parser.Results.n_mini_categ_skip; 
p_adj_method = parser.Results.p_adj_method; 
covariate = parser.Results.covariate; 
verbose = parser.Results.verbose; 

%% input checking

% Sanity check - there should be no nans. 
assert(all(~isnan(rdm(:))), ...
    'find_test_best_categ_model:nansInRDM', ...
    'There are NaNs in the RDM!'); 

% It can happen that the data and model matrix could be incompatible, in
% the sense that one is a similarity matrix (e.g. estimated using
% correlation across conditions), while the other is a dissimmilarity
% matrix (e.g. estimated using L2 norm across conditions). Since this
% funciton tests models that are based on similarity, the RDM of the data
% should have ones on the diagonal. 
assert(all(diag(rdm) == 1), ...
    'find_test_best_categ_model:dataNotSimilarityMatrix', ...
    'The input seems to be a dissimilarity matrix, but you need to provide a similarity matrix!'); 

%% prepare RDM

n_cond = size(rdm, 1); 

% get lower-triangular mask 
mask = tril(ones(n_cond, n_cond, 'logical')) & ~diag(ones(1, n_cond, 'logical')) ; 
mask_idx = find(mask); 

% set everything outside of lower triangular to NaN (just in case) 
rdm(~mask) = nan; 

%% get observed values

% generate all possible categorical models 
[models, all_models_info] = generate_categ_models(n_cond, n_categ, ...
                                'n_mini_categ_skip', n_mini_categ_skip); 

                            
% test all the models
[r, p, p_adj] = test_models(rdm, models, ...
                            'covariate', covariate, ...
                            'verbose', false);
                        
% add the returned values as columns to the table of models                            
all_models_info.r = r; 
all_models_info.p = p; 
all_models_info.p_adj = p_adj; 

% order the models from best to worst 
[~, idx] = sort(all_models_info.r, 'descend');
all_models_info = all_models_info(idx, :);
models = models(idx);
     
% get info about the best model 
r_obs = all_models_info{1, 'r'}; 
best_model = models{1};
best_model_info = all_models_info(1, :); 

%% do permutations to test significance of the best model

% compute null distribution of corelation coefficients 
r_null_dist = nan(1, n_perm); 
p_perm = nan; 

if n_perm > 0

    parfor i_perm=1:n_perm

        if mod(i_perm, 100) == 0 && verbose
            fprintf('perm %d/%d\n', i_perm, n_perm); 
        end

        % randomly shuffle the values on lower triangular of RDM1 
        shuffled_idx = randsample(mask_idx, length(mask_idx), false); 
        shuffled_rdm = rdm; 
        shuffled_rdm(mask_idx) = rdm(shuffled_idx); 

        % compute the correlation 
        r_perm = test_models(shuffled_rdm, ...
                        models, ...
                        'covariate', covariate, ...
                        'verbose', false); 

        % keep the highest correlation value
        r_null_dist(i_perm) = max(r_perm); 

    end

    % get pvalue from the null distribtuion 
    % see: https://stats.stackexchange.com/a/112456
    B = sum(r_null_dist > r_obs); 
    M = n_perm; 
    p_perm = (B+1) / (M+1); 

    if strcmpi(p_adj_method, 'bonf')
        p_perm = min(p_perm * length(models), 1); 
    end

end

