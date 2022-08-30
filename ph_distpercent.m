function [ ph_percent ] = ph_distpercent( ph_vowl, ph_cons, options )


nvowl = size(ph_vowl, 1);
ncons = size(ph_cons, 1);
nphoneme = nvowl + ncons;
ph_percent = zeros(nphoneme);
ph_vowl_sig = ph_vowl(:,2:(end-1));
ph_cons_sig = ph_cons(:,2:(end-1));
vowl_percent = ((nvowl:-1:1)/nvowl)*100;
cons_percent = ((ncons:-1:1)/ncons)*100;
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
for i = 1:ncons
   ph_code = ph_cons_code(i);       
   dist = zeros(ncons,1);
   for j = 1:ncons
      dist(j) = norm((ph_cons_sig(i,:)-ph_cons_sig(j,:)).*options.formant_weightscons'); 
   end
   [~,isort] = sort(dist);
   ph_code_sort = ph_cons_code(isort);
   for j = 1:ncons
       ph_percent(ph_code,ph_code_sort(j)) = cons_percent(j);
   end
end
end

