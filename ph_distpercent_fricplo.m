function [ ph_percent ] = ph_distpercent( ph_vowl, ph_cons, options, cons_fricplo_rows )


nvowl = size(ph_vowl, 1);
ncons = size(ph_cons, 1);
ncons_fricplo = length(cons_fricplo_rows);
cons_nonfricplo_rows = setdiff(1:ncons, cons_fricplo_rows);
ncons_nonfricplo = length(cons_nonfricplo_rows);
nphoneme = nvowl + ncons;
ph_percent = zeros(nphoneme);
ph_vowl_sig = ph_vowl(:,2:(end-1));
ph_cons_sig = ph_cons(:,2:(end-1));
vowl_percent = ((nvowl:-1:1)/nvowl)*100;
cons_fricplo_percent = ((ncons_fricplo:-1:1)/ncons_fricplo)*100;
cons_nonfricplo_percent = ((ncons_nonfricplo:-1:1)/ncons_nonfricplo)*100;


% calculate probabilities for vowels
for i = 1:nvowl
   ph_code = ph_vowl(i,1);
   dist = zeros(nvowl,1);
   for j = 1:nvowl
      dist(j) = norm((ph_vowl_sig(i,:)-ph_vowl_sig(j,:)).*options.formant_weightsvowl'); 
   end
   [~,isort] = sort(dist);
   ph_code_sort = ph_vowl(isort,1);
   for j = 1:nvowl
       ph_percent(ph_code,ph_code_sort(j)) = vowl_percent(j);
   end
end
% calculate probabilities for consonants
ph_cons_code = ph_cons(:,1);
ph_cons_code(ph_cons_code == -8) = 45;

ph_fricplo_code = ph_cons_code(cons_fricplo_rows);
ph_fricplo_sig  = ph_cons_sig(cons_fricplo_rows,:);
for i = 1:ncons_fricplo
   ph_code = ph_fricplo_code(i);       
   dist = zeros(ncons_fricplo,1);
   for j = 1:ncons_fricplo
      dist(j) = norm((ph_fricplo_sig(i,:)-ph_fricplo_sig(j,:)).*options.formant_weightscons'); 
   end
   [~,isort] = sort(dist);
   ph_code_sort = ph_fricplo_code(isort);
   for j = 1:ncons_fricplo
       ph_percent(ph_code,ph_code_sort(j)) = cons_fricplo_percent(j);
   end
end

ph_nonfricplo_code = ph_cons_code(cons_nonfricplo_rows);
ph_nonfricplo_sig  = ph_cons_sig(cons_nonfricplo_rows,:);
for i = 1:ncons_nonfricplo
   ph_code = ph_nonfricplo_code(i);       
   dist = zeros(ncons_nonfricplo,1);
   for j = 1:ncons_nonfricplo
      dist(j) = norm((ph_nonfricplo_sig(i,:)-ph_nonfricplo_sig(j,:)).*options.formant_weightscons'); 
   end
   [~,isort] = sort(dist);
   ph_code_sort = ph_nonfricplo_code(isort);
   for j = 1:ncons_nonfricplo
       ph_percent(ph_code,ph_code_sort(j)) = cons_nonfricplo_percent(j);
   end
end
end

