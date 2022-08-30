function [ FricPlo ] = FricPlo_check( ints1_ph, ints2_ph )


ll = length(ints1_ph);
l1 = floor(ll/4); % use only center half portion, 
l2 = ll-l1;
if( mean(ints2_ph(l1:l2)) >= .96*mean(ints1_ph(l1:l2)))
    FricPlo = 1;
else
    FricPlo = 0;
end

end

