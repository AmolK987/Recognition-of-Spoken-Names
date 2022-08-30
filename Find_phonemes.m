function Phoneme_list = Find_phonemes(id,phoneme2spel)
%
% function Phoneme_list = Find_phonemes(id,phoneme2spel)
%   Finds phonemes in the phoneme with spelling chart that correspond to one spelling
%
% Inputs
%   id - One spelling id
%   phoneme2spel - chart of phonemes with spelling
%
% Outputs
%   Phoneme_list - List of phonemes for one spelling id
%

[Phoneme_list,~] = find(abs(phoneme2spel - id)< 0.1);


end

