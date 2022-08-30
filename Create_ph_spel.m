function [ phoneme2spel ] = Create_ph_spel( ph_file, max_col )
%
% function [ phoneme2spel ] = Create_ph_spel( ph_file, max_col )
%   create the phoneme - spelling table from file
%
% Inputs
%   ph_file - name of file with phoneme table - string
%   max_col - Maximum number of spellings for any phoneme
%
% Outputs
%   phoneme2spel - matrix with phonemes in row and corresponding spelling in column
%

phon_num = 44; % number of phonemes

phoneme2spel = - ones(phon_num, max_col);

fid = fopen(ph_file);
for i = 1:phon_num
    ph_line = fgetl(fid);
    if ischar(ph_line)
        [ph_id,ph_sp_str] = fileline2phcode(ph_line);
        for n = 1:length(ph_sp_str)
            phoneme2spel(ph_id, n) = str2id(ph_sp_str{n}); 
        end
    end
end
fclose(fid);
end

