function [ formant_int, solid_formant_frac, formant_freq ] = phoneme_feat( sdb_norm, formants, f_band, f, solid_phoneme, mean_phoneme )
%
% function [ formant_int, solid_formant_frac ] = phoneme_feat( sdb_norm, formants, f_band, f, solid_phoneme, mean_phoneme )
%   create bands for frequency ranges of plosive "explosions"
%
% Inputs
%   sdb_norm - normalized spectrogram for the current phoneme
%   formants - the formants from Praat for the current phoneme
%   f_band - width of bands around each found formant
%   f - frequency vector
%   solid_phoneme - matrix of formant frequencies for when there is no noise
%
% Outputs
%   formant_int - vector with the intensities of each formant in the current phoneme
%   solid_formant_frac - 1 - fraction of phoneme duration with noise in each formant
%   formant_freq - vector with frequencies of all four formants at the midpoint of the phoneme
%

formant_int = zeros(size(formants,2),1);
solid_formant_frac = sum(solid_phoneme ~= 0, 1);
solid_formant_frac = solid_formant_frac'/size(solid_phoneme,1);
mid_row = round(size(mean_phoneme,1)/2);
formant_freq = mean_phoneme(mid_row,:)';

for i = 1:size(sdb_norm,2)
   f_curr = formants(i,:);
   for ii = 1:length(f_curr)
       F_top = f_curr(ii) + f_band;
       F_btm = f_curr(ii) - f_band;
       formant_int(ii) = formant_int(ii) + sum(sdb_norm((f >= F_btm) & (f<= F_top),i));
   end
end
formant_int = formant_int/size(sdb_norm,2); % getting average value with respect towards length of phoneme
end

