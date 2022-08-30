function [ match_percent ] = Match_addresbk( addresbk, found_name, ph_percent )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nnames = length(addresbk);
match_percent = zeros(nnames,1);
found_ph = length(found_name);
exlength_penalty = 10;


for i = 1:nnames
    name_match = addresbk{i};
    [ double_flag, trunc_names ] = double_cons( name_match );
    if ~double_flag
        tot_penalty = exlength_penalty*abs(length(name_match)-found_ph);
        nmatch = min([length(name_match) found_ph]);
        for j = 1:nmatch
            match_percent(i) = match_percent(i) + ph_percent(found_name(j),name_match(j));
        end
        match_percent(i) = match_percent(i)/nmatch - tot_penalty;
    else
        temp_match_percent = zeros(length(trunc_names),1);
        for l = 1:length(trunc_names)
            name_match = trunc_names{l};
            tot_penalty = exlength_penalty*abs(length(name_match)-found_ph);
            nmatch = min([length(name_match) found_ph]);
            for j = 1:nmatch
                temp_match_percent(l) = temp_match_percent(l) + ph_percent(found_name(j),name_match(j));
            end
            temp_match_percent(l) = temp_match_percent(l)/nmatch - tot_penalty;
            
        end
        match_percent(i) = max(temp_match_percent);
    end
end

end

