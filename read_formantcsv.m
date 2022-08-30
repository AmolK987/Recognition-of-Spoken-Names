function [t,formants,pitch,intensity] = read_formantcsv(fname,dt,intensity_thres,max_intensity)
% read formants csv file from .csv file


ff = fopen(fname,'r');

% fline = ' ';
nformants = 4;

max_rows = 2000;
t = zeros(max_rows,1);
formants = zeros(max_rows,nformants);
pitch = zeros(max_rows,1);
intensity = zeros(max_rows,1);

icount = 1;

fline = fgetl(ff); % skip header line

while ischar(fline)
    fline = fgetl(ff);
    
    fline = strrep(fline,'--undefined--',' 0 ');
    
    if ~ischar(fline)
        break
    end
    
%     if ~isempty(strfind(fline,'undefined'))
%         if icount == 1
%             t(icount) = 0;
%         else
%             t(icount) = t(icount-1)+dt;
%         end
%         formants(icount,:) = 0;
%     else
        tmp = sscanf(strrep(fline,',',' '),'%f');
        
        t(icount) = tmp(1);
        formants(icount,:) = tmp(1+(1:nformants))';
        pitch(icount) = tmp(nformants+2);
        intensity(icount) = tmp(nformants+3);
        
%         if intensity(icount) < intensity_thres
%             formants(icount,:) = 0;
%         end
%     end
    
    icount = icount+1;
end

if max_intensity > 0
    intensity_ratio = max_intensity / max(intensity);
    intensity = intensity * intensity_ratio;
end

formants(intensity<intensity_thres,:) = 0;


t = t(1:(icount-1));
formants = formants(1:(icount-1),:);
pitch = pitch(1:(icount-1),:);
intensity = intensity(1:(icount-1),:);

fclose(ff);
