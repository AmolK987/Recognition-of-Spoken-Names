function [ spel2phoneme,spel_list ] = Create_spel_ph( phoneme2spel, max_col)
%
% function [ spel2phoneme,spel_list ] = Create_spel_ph( phoneme2spel, max_col)
%   create the Spelling - Phoneme table from chart of phonemes with spellings
%
% Inputs
%   phoneme2spel - chart of phonemes with corresponding spellings (see Create_ph_spel)
%   max_col - Maximum number of phonemes for any spelling
%
% Outputs
%   spel2phoneme - matrix with spelling(id) in row and corresponding phonemes in column
%   spel_list - cell array of unique spelling strings
%

temp = phoneme2spel(:); % create 1 vector with all columns
temp2 = temp((temp>0)); % discard any -1, keep only spelling ids
temp2u = unique(temp2);   % identify only unique spelling ids and in ascending order

spel_num = length(temp2u);

spel_list = cell(spel_num,1);

spel2phoneme = -ones(spel_num, max_col);

spel2phoneme(:,1) = temp2u;

for i = 1:spel_num
   id = spel2phoneme(i,1);
   
   Phoneme_list = Find_phonemes(id,phoneme2spel);
   
   spel2phoneme(i,(2:(length(Phoneme_list)+1))) = Phoneme_list;
   
   spel_list{i} = id2str(id);
%    for p = 1:length(Phoneme_list)
%       spel2phoneme(i,p+1) = Phoneme_list(p); 
%    end
end


end

