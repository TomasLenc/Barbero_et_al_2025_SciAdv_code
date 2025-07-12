function [models, tbl] = generate_categ_models(n_cond, n_categ, varargin)
% This function generates all possible categorical models for the given
% number of conditions and categories. 
% 
% Parameters
% ----------
% n_cond : int
%     Number of unique conditions in the tested slice through the rhythm
%     space. 
% n_categ : int {2, 3, [2, 3]}
%     Number of categories to fit in the model. Either 2, 3, or [2, 3] will
%     generate all possible models that have 2 AND 3 categories.
% n_mini_categ_skip : int, optional, default=0
%     Models that include categories with <= `n_mini_categ_skip` pixels
%     will be excluded.
% 
% Returns 
% -------
% models : cell array 
%     Each element is a `n_cond` x `n_cond` matrix corresponding to one
%     categorical model. Observed correlation coefficient.
% tbl : table 
%     Matlab table where each row corresponds to one model, and describes
%     its properties. 
% 

parser = inputParser; 

addParameter(parser, 'n_mini_categ_skip', 0); 

parse(parser, varargin{:}); 

n_mini_categ_skip = parser.Results.n_mini_categ_skip; 

%%

col_names = {'n_categ', 'bound1', 'bound2'}; 
tbl = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 

models = {}; 

categ_bounds = [n_mini_categ_skip+1 : (n_cond - n_mini_categ_skip - 1)]; 


%% get 2-category models 

if ismember(2, n_categ)
    
    for i=1:length(categ_bounds)
        
        new_row = [{2}, {categ_bounds(i)}, {nan}]; 
        tbl = [tbl; new_row]; 

        models{end+1, 1} = get_2categ_model(n_cond, categ_bounds(i)); 

    end

end


%% get 3-category models

if ismember(3, n_categ)

    bounds1 = categ_bounds(1 : end - 1 - n_mini_categ_skip); 
    
    for i1=1:length(bounds1)

        bounds2 = [bounds1(i1) + n_mini_categ_skip + 1 : ...
                   (n_cond - n_mini_categ_skip - 1)]; 

        for i2=1:length(bounds2)

            models{end+1, 1} = get_3categ_model(n_cond, bounds1(i1), bounds2(i2)); 

            new_row = cellfun(@num2cell, {3, bounds1(i1), bounds2(i2)}); 
            tbl = [tbl; new_row]; 
            
        end

    end

end
