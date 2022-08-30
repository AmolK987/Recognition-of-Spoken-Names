function [ str_id ] = str2id( strg )
%
%function [ str_id ] = str2id( strg )
%   Take any string of letters a-z and convert to a unique number
%
% Inputs
%   strg - string of letters
%
% Outputs
%   str_id - unique number for string
%

strg    = lower(strg); % ensure all characters in lower case
str_len = length(strg);

str_id = 0;

for n = 1:str_len
   ch = strg(n);
   str_id = str_id + (uint32(char2id(ch)) * 27^(str_len - n));
   
end


end

