function [ First_pseudo, Last_pseudo, RealF_pseudo, RealL_pseudo ] = create_pseudo( First_str, Last_str, Real_first, Real_last )
%
% function [ First_pseudo, Last_pseudo, RealF_pseudo, RealL_pseudo ] = create_pseudo( First_str, Last_str, Real_first, Real_last )
%   Take First and Last names (real and produced) and convert into pseudo-spellings
% 
%   Inputs
%    First_str - String for the produced first name from phonemes
%    Last_str - String for the produced last name from phonemes
%    Real_first - True spelling of first name from dictionary - string
%    Real_last - True spelling of last name from dictionary - string
%
%   Outputs
%    First_pseudo - Pseudo-spelling of produced first name
%    Last_pseudo - Pseudo-spelling of produced last name
%    RealF_pseudo - Pseudo-spelling for the real first name
%    RealL_pseudo - Pseudo-spelling for the real last name
%
ph_chart = Create_ph_spel('phoneme_chart.csv', 7);

phf = spelconvph(First_str);
phl = spelconvph(Last_str);
phRf = spelconvph(Real_first);
phRl = spelconvph(Real_last);

First_pseudo = ph2pseudo(phf, ph_chart);
Last_pseudo = ph2pseudo(phl, ph_chart);
RealF_pseudo = ph2pseudo(phRf, ph_chart);
RealL_pseudo = ph2pseudo(phRl, ph_chart);
end

