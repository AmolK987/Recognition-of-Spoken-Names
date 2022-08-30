function [ compf ] = Compare_pseudo( f_str, real_f )
%
% function [ compf ] = Compare_pseudo( f_str, real_f )
%   Compare the pseudo-spelling to the real spelling from dictionary
%
% Inputs
%   f_str - pseudo-spelling string
%   real_f - pseudo spelling from dictionary
%
% Outputs
%   compf - percent for match between the pseudo-spelling and pseudo-spelling from dictionary
%

compf = 100;
lf = length(f_str);
lrf = length(real_f);
sl = min(lf,lrf);
ml = max(lf,lrf);

wrong =  abs( lf - lrf);
for ii = 1:sl
    lf_ii = f_str(ii);
    lrf_ii = real_f(ii);
    if (lrf_ii ~= lf_ii)
        wrong = wrong + 1;
    end
end


compf = compf - (wrong/ml)*100;


end

