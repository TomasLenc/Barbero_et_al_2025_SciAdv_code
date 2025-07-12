function model=get_3categ_model(n, bound1, bound2)
% prepare categorical model RDM with 3 categories based on requested location 
% of the boundaries

assert(length(n) == 1, 'The number of requested conditions must be a single integer'); 
assert(length(bound1) == 1, 'The first boundary must be a single integer'); 
assert(length(bound2) == 1, 'The second boundary must be a single integer'); 

assert(bound1 >= 1 && bound1 < n, ...
      'get_3categ_model:firstBoundOutOfRange', ...
      'The first boundary must be between 1 and N_conditions-1'); 
  
assert(bound2 >= 1 && bound2 < n, ...
      'get_3categ_model:secondBoundOutOfRange', ...
      'The first boundary must be between 1 and N_conditions-1'); 
  
assert(bound2 > bound1, ...
      'get_3categ_model:secondBoundSmallerThanFirst', ...
      'The second boundary must be larger than the first one'); 
  

model = zeros(n, n); 
for i_row=1:n
    for i_col=1:n
        
        if i_row <= bound1 && i_col <= bound1
            model(i_row, i_col) = 1; 
        end
        
        if i_row > bound1 && i_col > bound1 && ...
                i_row <= bound2 && i_col <= bound2
            model(i_row, i_col) = 1; 
        end
        
        if i_row > bound2 && i_col > bound2
            model(i_row, i_col) = 1; 
        end
    end
end


