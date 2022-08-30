function [ char_id ] = char2id( charac )
%
% function [ char_id ] = char2id( charac )
%   convert a character a-z to number code 1-26
%
% Inputs
%  char : input character, a-z
%
% Outputs
%  char_id: corresponding number code a=1, z = 26

char_id = uint8(charac) - 96;

end

