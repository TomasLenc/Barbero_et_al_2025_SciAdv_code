function model=get_2categ_model(n, bound)
% prepare categorical model RDM based on requested boundary location 

assert(bound >= 1 && bound < n, ...
      'get_2categ_model:boundOutOfRange', ...
      'The first boundary must be between 1 and N_conditions-1'); 
  
model = zeros(n, n); 
for i_row=1:n
    for i_col=1:n
        
        if i_row <= bound && i_col <= bound
            model(i_row, i_col) = 1; 
        end
        if i_row > bound && i_col > bound
            model(i_row, i_col) = 1; 
        end
    end
end