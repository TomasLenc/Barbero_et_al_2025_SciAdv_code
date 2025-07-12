function [partials, dp_quad, dp_cub, all_f] = get_partials(trial_dur, partials)
% This function finds the closest possible frequencies for each requested
% partial that complete an exact number of cycles within the duration of the
% trial. It will also print a report to the command line. 
%
%

% partials = [234 311 577];
% trial_dur = 30 * 0.900; 

% find closest frequencies that will have integer N of cycles in the trial AND
% will be an integer ratio with the sampling rate 
for i_f=1:length(partials)

    % find closest integer ratios to trial duration 
    target_period = 1 / partials(i_f); 

    n_cycles = trial_dur / target_period; 
    
    closest_f = round(n_cycles) / trial_dur; 
    
    partials(i_f) = closest_f; 

end

% dp_quad = [partials(2)-partials(1), ...
%            partials(3)-partials(2),...
%            partials(3)-partials(1)]; 
    
A = partials - partials'; 
m  = triu(true(size(A)), 1);
dp_quad = sort(A(m)'); 
       
if length(partials) == 2
    dp_cub = 2*partials(1)-partials(2);
elseif length(partials) == 3
    dp_cub = [2*partials(1)-partials(2), ...
              2*partials(2)-partials(3)]; 
else 
    error('there must be 2 or 3 partials - you asked for %d', length(partials)); 
end

% if we have some partials in a harmonic ratio, the DPs calculated above
% will overlap with the partials... let's make sure this doesn't happen
dp_quad = dp_quad(dp_quad > 0); 
dp_quad = dp_quad(~ismember(dp_quad, partials)); 

dp_cub = dp_cub(dp_cub > 0); 
dp_cub = dp_cub(~ismember(dp_cub, partials)); 


%%

fprintf('\n\n=======================================================================\n'); 

fprintf('\npartials:   %s \n', vec_to_str(partials, 'sep', '   ', 'format', '%.3f')); 

if length(partials) == 2
    
    fprintf('\npartial ratios: \n2/1 = \t%.3f\n\n', ...
            partials(2) / partials(1) ...
            ); 
        
    
elseif length(partials) == 3
    
    fprintf('\npartial ratios: \n2/1 = \t%.3f \n3/2 = \t%.3f \n3/1 = \t%.3f\n\n', ...
            partials(2) / partials(1), ...
            partials(3) / partials(2), ...
            partials(3) / partials(1) ...
            ); 
        
        

end

fprintf('quadratic DPs (f2-f1): \t%.3f\n', ...
    vec_to_str(dp_quad, 'sep', '\t', 'format', '%.3f')); 

fprintf('cubic DPs (2*f1-f2): \t%.3f\n', ...
    vec_to_str(dp_cub, 'sep', '\t', 'format', '%.3f')); 
  
    
all_f = sort([partials, dp_quad, dp_cub]'); 
all_dp = sort([dp_quad, dp_cub]'); 

fprintf('\nall DPs sorted: \n\n'); 
disp(all_dp)

fprintf('\nall frequencies sorted: \n\n'); 
disp(all_f)

fprintf('\n[-5Hz +5Hz] range around each frequency: \n\n'); 
disp([all_f-5, all_f+5])

fprintf('\nis multiple of 50Hz included within the +- 5Hz interval?: \n\n'); 
disp([mod((all_f-5),50) > mod(all_f+5,50)])

fprintf('\nmodulus trial_duration/partial_period  (sanity check) \n'); 
for i=1:length(partials)
    fprintf('%.2f Hz : %.11f\n', partials(i), mod(trial_dur / (1/(partials(i))), 1)); 
end
    
fprintf('\nmodulus trial_duration/crossmod_period (sanity check) \n'); 
for i=1:length(all_dp)
    fprintf('%.2f Hz : %.11f\n', all_dp(i), mod(trial_dur / (1/(all_dp(i))), 1)); 
end
    
    
    
    
    
    
    
    
    