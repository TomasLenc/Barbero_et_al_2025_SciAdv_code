function labs = get_ratio_labs(ratios, varargin)

parser = inputParser; 

addParameter(parser, 'n_decim', 2); 
addParameter(parser, 'no_zero', false); 

parse(parser, varargin{:}); 

n_decim = parser.Results.n_decim; 
no_zero = parser.Results.no_zero; 


if size(ratios, 2) == 1
  labs = cellfun(@(x) num2str(x, sprintf('%%.%df', n_decim)), num2cell(ratios), 'uni', 0);     
else
  labs = cellfun(@(x) vec_to_str(x ./ min(x), 'format', sprintf('%%.%dg', n_decim), 'sep', ', '), ...
                 num2cell(ratios, 2), 'uni', 0); 
end

if no_zero
    labs = cellfun(@(x) regexprep(x, '^0.', '.'), labs, 'uni', 0); 
end