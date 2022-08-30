function spelstr = id2str(spelcode)
%
% function spelstr = id2str(spelcode)
%   convert a numeric value of string back to string - opposite of str2id
%  
% Inputs
%   spelcode : a numeric value of the spelling string
%  
% Outputs
%   spelstr  : string corresponding to spelcode
%

dspelcode = double(spelcode);

if dspelcode < 1.5 % special case when spelcode = 1 (for a)
    n = 1;
else
    n = ceil(log(dspelcode)/log(27));
end;

rem = dspelcode;

for ii = 1:n,
    
    if ii < n
        temp = floor(rem/(27^(n-ii))); % get the number for ii-th character
    else
        temp = rem;
    end;
    spelstr(ii) = char(temp+96);   % convert number to char
    
    if ii < n
        rem = rem - temp*27^(n-ii);
    end;
end

