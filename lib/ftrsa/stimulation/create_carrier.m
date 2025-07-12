function [s, t] = create_carrier(partials, duration, fs)

N = round(duration * fs); 

s = zeros(1, N); 

t = [0 : N-1] / fs; 

for i_f=1:length(partials)
    
    s = s + sin(2 * pi * partials(i_f) * t); 
    
end

s = s ./ max(abs(s)); 