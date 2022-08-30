function [ ph_trans, ph_start, locs1, locs2, ph_trans_ind ] = Calc_ph_trans( int_filt_word, word_t, peakfnd_start, peak_dist, peak_prominence, first_consonant, pitch_start )

[pks1,locs_ind1,width1,prom1]= findpeaks(int_filt_word,'MinPeakDistance',peak_dist,'NPeaks',4,'MinPeakProminence',peak_prominence);% maximum peaks
[pks2,locs_ind2,width2,prom2]= findpeaks(-int_filt_word,'MinPeakDistance',peak_dist,'NPeaks',4,'MinPeakProminence',peak_prominence); % minimum peaks

try_count = 1;
while (isempty(locs_ind2) || isempty(locs_ind1)) && (try_count <=4);
    peak_prominence = 0.75*peak_prominence;
    disp('trying with reduced prominence');
    [pks1,locs_ind1,width1,prom1]= findpeaks(int_filt_word,'MinPeakDistance',peak_dist,'NPeaks',4,'MinPeakProminence',peak_prominence);% maximum peaks
    [pks2,locs_ind2,width2,prom2]= findpeaks(-int_filt_word,'MinPeakDistance',peak_dist,'NPeaks',4,'MinPeakProminence',peak_prominence); % minimum peaks
    try_count = try_count+1;
end

locs1 = word_t(locs_ind1);
locs2 = word_t(locs_ind2);
if (first_consonant == 1) && (word_t(pitch_start) > locs1(1)) && (locs1(1) < locs2(1))% found a maximum before pitch start - remove it
   locs1 = locs1(2:end); 
   locs_ind1 = locs_ind1(2:end);
end
if (first_consonant == 0) && (isempty(locs2))
   locs2 = word_t(end);
   locs_ind2 = length(word_t);
end
if (first_consonant == 1) && (isempty(locs2) || (locs2(1) > word_t(pitch_start)) && (locs1(1) < locs2(1))) % found no mimimum before pitch start - add one
   min_pad_ind = round(mean([1, locs_ind1(1)]));
   min_pad = word_t(min_pad_ind);
   locs2 = [min_pad; locs2];
   locs_ind2 = [min_pad_ind; locs_ind2];
end
max_ind = find(locs1 >= peakfnd_start,1);
min_ind = find(locs2 >= peakfnd_start,1);
diff_int = [0;diff(int_filt_word)];
ph_curr = 1;
if locs1(max_ind) < locs2(min_ind)
   ph_start = 1; % max first meaning vowel
   peak_dir = -1;
else
    ph_start = 0; % min first meaning consonant
    peak_dir = 1;
end
pks_done = 0;
while ~pks_done
   if peak_dir == -1
%       [pks,locs]= findpeaks(-diff_int(locs_ind1(max_ind):locs_ind2(min_ind)),'NPeaks',1);
      [pks,locs]= min(diff_int(locs_ind1(max_ind):locs_ind2(min_ind)));
      temp_ind = locs_ind1(max_ind)+ (locs(1)-1);
      ph_trans(ph_curr) = word_t(temp_ind);
      ph_trans_ind(ph_curr) = temp_ind;
      max_ind = max_ind + 1;
      if max_ind > length(locs_ind1)
         pks_done = 1; 
      end
      peak_dir = 1; 
   else
%       [pks,locs]= findpeaks(diff_int(locs_ind2(min_ind):locs_ind1(max_ind)),'NPeaks',1);
      [pks,locs]= max(diff_int(locs_ind2(min_ind):locs_ind1(max_ind)));
      temp_ind = locs_ind2(min_ind)+ (locs(1)-1);
      ph_trans(ph_curr) = word_t(temp_ind);
      ph_trans_ind(ph_curr) = temp_ind;
      min_ind = min_ind + 1;
      if min_ind > length(locs_ind2)
         pks_done = 1; 
      end
      peak_dir = -1;
   end
   ph_curr = ph_curr + 1;
end
ph_trans = [ph_trans(:); word_t(end)];
ph_trans_ind = [ph_trans_ind(:);length(word_t)];
if peakfnd_start > word_t(1) % include first consonant if pitch gap is larger than the threshold
   ph_trans = [peakfnd_start; ph_trans];
   ii = find(word_t >= peakfnd_start,1);
   ph_trans_ind = [ii; ph_trans_ind];
   ph_start = 0; 
end
end

