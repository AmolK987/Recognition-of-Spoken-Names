function [ norm_sdbj ] = Normalizer( sdbj,max_sdb,broadbnd,f,ic,j1,j2,max_ang )
%
% function [ norm ] = Normalizer( sdbj,max_sdb,broadbnd,f )
%   Takes sdbj and finds points where data value decreases to locate and nomalize peaks
%
% Inputs
%   sdbj - each column of data from sdb
%   max_sdb - max value in sdb to return normalized value between 0-1
%   broadbnd - using broad band(1) or narrow band(0)
%   f - frequency vector
%   ic - index of column of sdb
%   j1 - index of first column included in band around ic
%   j2 - index of last column included in band around ic
%   max_ang - max angle to search over (in radians)
%
% Outputs
%   norm_sdbj - normalized matrix of 0 and 1
%

% filter signal first to remove noise

% if ~broadbnd
%    [b,a] = butter(2,0.4);
% else
%     [b,a] = butter(2,0.99);%0.2
% end
% alpha = 0.01;
% b = [1-alpha 0];
% a = [1 -alpha];
% 
% sdbjf1 = filtfilt(b,a,sdbj);
sdbjf1 = sdbj;

nang = 15; % odd no. of angles to search over, middle one is 0 deg
ang = linspace(-max_ang,max_ang,nang);

nrows = size(sdbj,1);
sdbjf = zeros(nrows,1);
for row = 1:nrows
    maxval = -10000;
    for ia = 1:nang
        angi = ang(ia);
        tangi= tan(angi);
        sum = 0;
        for i1 = ic:j2
            irow = round(row - (ic-i1)*tangi);
            if (irow > nrows) || (irow < 1)
                break
            end
            sum = sum+sdbjf1(irow,i1-j1+1);
        end
        if ic > 1
            for i1 = (ic-1):-1:j1
                irow = round(row - (ic-i1)*tangi);
                if (irow > nrows) || (irow < 1)
                    break
                end
                sum = sum+sdbjf1(irow,i1-j1+1);
            end
        end
        if sum > maxval
            maxval = sum;
        end
    end
    sdbjf(row) = maxval;
end
    
%sdbjf  = sum(sdbjf1,2);

MPDlo = 100; % frequency separation of peaks at low freq
MPDhi = 300; % freq separation of peaks at high freq
% MPDlo = 25; % frequency separation of peaks at low freq
% MPDhi = 150; % freq separation of peaks at high freq
freq_thresh = 2000; % location(Hz) of separation between low and high frequencies
freq_row = find(f >= freq_thresh, 1);
if isempty(freq_row),
    freq_row = length(sdbjf)+1;
end
npeaks_lo = 3;
npeaks_hi = 3;
MPP_lo = 0.05;
MPP_hi = 0.25;
[~,pkf_lo] = findpeaks(sdbjf(2:freq_row-1),f(2:freq_row-1),'MinPeakDistance',MPDlo,'SortStr','descend','NPeaks',npeaks_lo,'MinPeakProminence',MPP_lo);
[~,pkf_hi] = findpeaks(sdbjf(freq_row:end),f(freq_row:end),'MinPeakDistance',MPDhi,'SortStr','descend','NPeaks',npeaks_hi,'MinPeakProminence',MPP_hi);
% [~,pkf_lo] = findpeaks(sdbjf(2:freq_row-1),f(2:freq_row-1),'MinPeakDistance',MPDlo,'NPeaks',npeaks_lo,'MinPeakProminence',MPP_lo);
% [~,pkf_hi] = findpeaks(sdbjf(freq_row:end),f(freq_row:end),'MinPeakDistance',MPDhi,'NPeaks',npeaks_hi,'MinPeakProminence',MPP_hi);
norm_sdbj = zeros(size(sdbjf));
pkf = [pkf_lo ; pkf_hi];

% figure(3); clf
% plot(sdbj); hold on; plot(sdbjf,'r--');
% disp('plotted');

% fsdbj = [0; diff(sdbjf)];
% 
% fsdbj_sh = [0; fsdbj(1:end-1)];
% 
% ii = find((fsdbj < 0) & (fsdbj_sh >= 0));
% 
if ~isempty(pkf)
    [comm,jj,ii]= intersect(pkf,f);
    norm_sdbj(ii-1) = sdbj(ii-1);
end

if ~broadbnd
    norm_sdbj(norm_sdbj < 0.65) = 0;
else
    thres_lo = 0.7;
    thres_hi = 0.1;
    thres    = [repmat(thres_lo,freq_row,1); repmat(thres_hi,nrows-freq_row,1)];
    norm_sdbj(norm_sdbj < thres) = 0; %0.4
end
