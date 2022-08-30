function [ pseudo_spel ] = ph2pseudo( ph_str, ph_chart )
%
%   function [ pseudo_spel ] = ph2pseudo( ph_str, ph_chart )
%       Takes string of phonemes(ph_str) and converts to pseudo-spelling(pseudo_spel)
%
% Inputs
%  ph_str - array of phonemes for corresponding name
%  ph_chart - chart of phonemes with corresponding spelling ids
%
% Outputs
%  pseudo_spel - string of Pseudo-spelling - only one spelling for each corresponding phoneme
%
x = length(ph_str);
pseudo_spel = '';

for ii = 1:x
   ph_id = ph_str(ii); 
   spel_id = ph_chart(ph_id, 1);
   spel = id2str(spel_id);
   pseudo_spel = [pseudo_spel spel]; 
end

end

