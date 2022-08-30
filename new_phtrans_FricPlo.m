function [ ph_trans, force_frstcons, ph_fricplo ] = new_phtrans_FricPlo( ph_trans, ph_trans_ind, ints1, ints2, first_consonant )

force_frstcons = 0;
remove_trans = zeros(length(ph_trans),1);
cons_sequence = zeros(length(ph_trans),1);
if(first_consonant ==1)
    cons_sequence(1:2:end) = 1;
else
    cons_sequence(2:2:end) = 1;
end
for i = 1:length(ph_trans)
    if cons_sequence(i) == 0
        if ( i == 1)
            [ FricPlo ] = FricPlo_check( ints1(1:ph_trans_ind(1)), ints2(1:ph_trans_ind(1)));
            if FricPlo == 1
                remove_trans(1) = 1;
                force_frstcons = 1;
                first_consonant= 1;
            end
            
        elseif ( i == length(ph_trans))
            [ FricPlo ] = FricPlo_check( ints1(ph_trans_ind(end-1):ph_trans_ind(end)), ints2(ph_trans_ind(end-1):ph_trans_ind(end)));
            if FricPlo == 1
                remove_trans(end-1) = 1;
            end
        else
            [ FricPlo ] = FricPlo_check( ints1(ph_trans_ind(i-1):ph_trans_ind(i)), ints2(ph_trans_ind(i-1):ph_trans_ind(i)));
            if FricPlo == 1
                remove_trans(i-1) = 1;
                remove_trans(i) = 1;
            end
        end
    end
end
ph_trans_orig = ph_trans;
ph_trans = ph_trans_orig(remove_trans==0); % keep the ones not marked to be removed
ph_trans_remove = ph_trans_orig(remove_trans==1); % ones marked to be removed
ph_fricplo = zeros(length(ph_trans),1);
if first_consonant
    ph_fricplo(1) = FricPlo_check( ints1(1:ph_trans_ind(1)), ints2(1:ph_trans_ind(1)));
end

cons_curr = ~first_consonant;
for i=2:length(ph_trans)
    if cons_curr
        i1 = ph_trans(i-1);
        i2 = ph_trans(i);
        if ~isempty(find((ph_trans_remove > i1) & (ph_trans_remove < i2),1))
            ph_fricplo(i) = 1;
        end
    end
    cons_curr = ~cons_curr;
end

