function p = my_gaus_mix(theta, x, varargin)
% theta:  mu, sigma, w
parser = inputParser; 
addParameter(parser, 'smooth', false); 

parse(parser, varargin{:}); 

do_smoothing = parser.Results.smooth; 

mu1 = theta(1); 
mu2 = theta(2); 
s1 = theta(3); 
s2 = theta(4); 
w1 = theta(5); 
w2 = theta(5); 

p = zeros(1, length(x)); 
    
p = p + w1 * 1/sqrt(2*pi*s1^2) * exp(-1/2 * ((x-mu1)/s1).^2); 
p = p + w2 * 1/sqrt(2*pi*s2^2) * exp(-1/2 * ((x-mu2)/s2).^2); 
    
p = p ./ sum(p); 

if do_smoothing
    p = p + 1e-5; 
    p = p ./ sum(p); 
end