function [ match_percent, name_ind ] = best_match_addresbk( found_name, addresbk, ph_percent )


% match for first names
[ first_percent ] = Match_addresbk( addresbk.first, found_name.first, ph_percent );

% match for last names
[ last_percent ] = Match_addresbk( addresbk.last, found_name.last, ph_percent );

match_percent.first = first_percent;
match_percent.last = last_percent;
match_percent.both = mean([first_percent last_percent],2);
[~,name_ind] = max(match_percent.both);
end

