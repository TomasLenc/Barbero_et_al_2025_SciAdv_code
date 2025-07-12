function s = insert_sound_events(s, idx, event, varargin)

parser = inputParser; 

addParameter(parser, 'vol', []); 
addParameter(parser, 'add', false); 

parse(parser, varargin{:}); 

vol = parser.Results.vol; 
add_sound = parser.Results.add; 

if ~isempty(vol) && length(vol) ~= length(idx)
   error('%d volume values but %d events', length(vol), length(idx));   
end

%%

for i=1:length(idx)
    
    if isempty(vol)
        gain = 1; 
    else
        gain = vol(i); 
    end
    
    if add_sound
        s(idx(i)+1 : idx(i)+length(event)) = ...
                        s(idx(i)+1 : idx(i)+length(event)) + gain * event; 
    else
        s(idx(i)+1 : idx(i)+length(event)) = gain * event; 
    end
    
end