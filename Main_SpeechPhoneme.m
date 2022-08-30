
clear all
warning off


% wav_fname = 'ssuresh_new.wav'; %71.7 yes
% wav_fname = 'amuppidi_new.wSheeav'; % 34.4 no, 72.3 yes
% wav_fname = 'glvov_new.wav'; % 49.4 no, 67.0 yes
% wav_fname = 'arvindkr_new.wav'; % 83.1 yes
% wav_fname = 'kyamaguchi_new.wav'; % 63.0 yes
% wav_fname = 'amolkr_new.wav'; % 58.9 yes
% wav_fname = 'ssuresh_new_papa.wav'; % 64.3 yes
% wav_fname = 'amuppidi_new_papa.wav'; % 28.1 no
% wav_fname = 'kyamaguchi_new_papa.wav'; % 61.1 yes
% wav_fname = 'arvindkr_new_papa.wav'; % 49.6 no
% wav_fname = 'glvov_new_papa.wav'; % 45.5 yes
% wav_fname = 'amolkr_new_papa.wav'; % 48.6 yes
% adoyle_new.wav - 18.5 yes
% bbrula_new.wav - 55.4 yes
% esarecha_new.wav - no - double vowel so only found mins
% jklamka_new.wav - 50.7 no
% jxue_new.wav - 38.9 no
% ncao_new.wav - 66.1 yes
% smanga_new.wav - 45.0 no

new_wav = dir('*_new.wav');
new_mom_wav = dir('*_new_mom.wav');
new_papa_wav = dir('*_new_papa.wav');
new_all_wav = [new_wav; new_papa_wav; new_mom_wav];

for i=1:length(new_all_wav)
    disp([int2str(i) '.   ' new_all_wav(i).name]);
end
wav_file_use = input('Enter no. of file to use from list : ');
if wav_file_use == 0
    wav_fnames = cell(length(new_all_wav),1);
    for i=1:length(new_all_wav)
        wav_fnames{i} = new_all_wav(i).name;
    end
else
    wav_fnames{1} = new_all_wav(wav_file_use).name;
    [data,fs] = audioread(wav_fnames{1});
    
    if(wav_file_use ~= 0)
        disp('-----------------------');
        disp(['Testing with file: ' wav_fnames{1}]);
        player = audioplayer(data, fs);
        figure(1); clf
        pause(1);
        t = (0:(1/fs):(length(data)-1)/fs);
        plot(t,data);
        grid on
        xlabel('Time(s)');
        title('Waveform');
        drawnow
        pause(1);
        %%play(player);
    end
end

best_match_percent = zeros(length(wav_fnames),3);
best_match_name = cell(length(wav_fnames),1);
correct_name = cell(length(wav_fnames),1);
first_cons = zeros(length(wav_fnames),2);


% weights over frequency for intensity calculation
int_weight_freq = [0 1; 4000 1; 8000 1; 8001 0; 100000 0];

error_flag = zeros(length(wav_fnames),1);

tic

for inames = 1:length(wav_fnames)
    
    try
        
    wav_fname = wav_fnames{inames};
    correct_name{inames} = wav_fname;
    
    % praat script parameters
    praat_tstart = 0; % start time for speech
    praat_tend = 6.9; % end time for speech
    praat_dt = 0.005; % sample time for praat results from script
    praat_intensity_thres = 54;%54; %55; %threshold for returning nonzero formants
    
    [data,fs] = audioread(wav_fname);
    data = data*4e4;
    
    timewin = 0.005;
    timeover = 0.5*timewin;
    
    win = round(timewin*fs);
    overlap = round(timeover*fs);
    % plotopt = 1;
    nfft = 512;
    [s1,sdb_f,sdb_t,p,fc,tc]= spectrogram(data,win,overlap,nfft,fs);
    int_wt_freq_sdbf = interp1(int_weight_freq(:,1),int_weight_freq(:,2),sdb_f);
    sdb_t = sdb_t(:);
    s1 = abs(s1);
    sdb = max(20*log10(s1),0);
    sdb_norm = sdb/max(max(max(sdb)),70);
    low_limit = 30; % anything below this will be white
    spec_cmap = ones(64,3);
    temp = linspace(1,0,(64-low_limit+1))';
    for j = 1:length(temp)
        spec_cmap(low_limit+(j-1), :) = temp(j);
    end
    
    delete('praat_formants.csv');
    pause(3);
    
    %!"Praat.exe" --run "testamolkr_formants.praat"
    % system('"Praat.exe" --run "testamolkr_formants.praat"');
    praat_timestr = [num2str(praat_tstart) ' ' num2str(praat_tend) ' ' num2str(praat_dt)];
    praat_fname = wav_fname;
    praat_cmd = ['"Praat.exe" --run calc_praat_formants.praat "' praat_fname '" ' praat_timestr];
    system(praat_cmd);
    pause(4);
    
    max_intensity = 80; % normalize peak intensity in recording to this value, if 0, then don't normalize
    
    [formant_t,formants,pitch,intensity]= read_formantcsv('praat_formants.csv',praat_dt,praat_intensity_thres, max_intensity);
    
    if formant_t(end) < sdb_t(end)
        last_row = find(sdb_t > formant_t(end),1)-1;
        sdb_t = sdb_t(1:last_row);
        sdb   = sdb(:,1:last_row);
        sdb_norm = sdb_norm(:,1:last_row);
    end
    formants = interp1(formant_t,formants,sdb_t);
    pitch = interp1(formant_t,pitch(:),sdb_t);
    intensity = interp1(formant_t,intensity(:),sdb_t);
    formant_t = sdb_t;
    intensityslope2 = [0; 0; diff(diff(intensity))]; % used to check first_consonant if pitch gap is too small
    intensityslope1 = [0; diff(intensity)]; % used to check first_consonant if pitch gap is too small
    
    
    med_window = 11;
    word_gapt_threshold = 0.11; % mimimum gap between words in seconds, below this gap - merge these words
    word_gapi_threshold = round(word_gapt_threshold/(formant_t(2)-formant_t(1)));
    word_lengtht_thresh = 0.3; % discard any word shorter than this
    word_lengthi_thresh = round(word_lengtht_thresh/(formant_t(2)-formant_t(1)));
    [ zstart, zend ] = Word_start_end( formants, med_window, word_gapi_threshold, word_lengthi_thresh );
    nword = length(zstart);
    zend = zend - 1;
    
    
    figure(1); clf
    mesh(sdb_t,sdb_f,sdb_norm); title('Spectrogram and Formants'); view([0 90]); colormap(spec_cmap)
    hold on;
    [x,y] = meshgrid(formant_t,sdb_f);
    ax =axis;
    ax(4) = 8000;
    axis(ax)
    for i=1:length(zstart)
        plot3([formant_t(zstart(i)); formant_t(zstart(i))],ax(3:4),[1.1; 1.1],'r--','linewidth',3.0);
        plot3([formant_t(zend(i)); formant_t(zend(i))],ax(3:4),[1.1; 1.1],'r--','linewidth',3.0);
    end
    xlabel('Time (s)'); ylabel('Freq (Hz)');
    %z = zeros(size(x));
    ncol = size(formants,2);
    linestdot = {'b.' 'g.' 'r.' 'm.'};
    linestsolid = {'b-' 'g-' 'r-' 'm-'};
    linestcircle = {'bo' 'go' 'ro' 'mo'};
    
    for i=1:ncol
        fvec = formants(:,i);
        zvec = zeros(size(fvec));
        zvec(fvec>0) = 1.1;
        plot3(formant_t,fvec,zvec,linestdot{i})
    end
    %plot(t,formants,'.');
    
    %% phoneme signatures
    
    read_Ph_excel = 0;
    if read_Ph_excel == 1
        xlfname = 'Broadband_ph_NewMic.xlsx';
        xlsheet = 'Phonemes';
        [ph_code,mF1,mF2,mF3,sF1,sF2,sF3,eF1,eF2,eF3,F1int,F2int,F3int,F1noise,F2noise,F3noise,pitchy]= read_phoneme_info(xlfname,xlsheet);
        save('ph_sig','ph_code','mF1','mF2','mF3','sF1','sF2','sF3','eF1','eF2','eF3','F1int','F2int','F3int','F1noise','F2noise','F3noise','pitchy');
    else
        load ph_sig
    end
    % phinfo = struct('ph_code', ph_code, 'mF1', mF1, 'mF2', mF2, 'mF3', mF3, 'sF1', sF1, 'sF2', sF2, 'sF3', sF3,...
    %     'eF1', eF1, 'eF2', eF2, 'eF3', eF3, 'F1int', F1int, 'F2int', F2int, 'F3int', F3int, 'F1noise', F1noise, 'F2noise', F2noise,...
    %     'F3noise', F3noise, 'pitchy', pitchy);
    consrows = 1:26;
    vowlrows = 27:45;
    temp = [ph_code mF1 mF2 mF3 sF1 sF2 sF3 eF1 eF2 eF3 F1int F2int F3int F1noise F2noise F3noise pitchy];
    ph_cons = temp(consrows,:);
    ph_cons_scale = mean(ph_cons(:,2:(end-1)), 1);
    ph_cons_scale((end-5):end) = 1;
    ph_cons(:,2:(end-1)) = ph_cons(:,2:(end-1))./repmat(ph_cons_scale,length(consrows),1);
    ph_cons_norm = zeros(length(consrows),1);
    for i = 1:length(consrows)
        ph_cons_norm(i) = norm(ph_cons(i, 2:(end-1)));
    end
    ph_vowl = temp(vowlrows,:);
    ph_vowl_scale = mean(ph_vowl(:,2:(end-1)), 1);
    ph_vowl_scale((end-5):end) = 1;
    ph_vowl(:,2:(end-1)) = ph_vowl(:,2:(end-1))./repmat(ph_vowl_scale,length(vowlrows),1);
    ph_vowl_norm = zeros(length(vowlrows),1);
    for i = 1:length(vowlrows)
        ph_vowl_norm(i) = norm(ph_vowl(i, 2:(end-1)));
    end
    
    use_praat_int_phtrans = 0;

    calcph_options.prev_ph_type = 0; % no previous phoneme
    %     calcph_options.noise_count_thresh = 1.5;  % number of first three formants that have noise must exceed this limit to mark a phoneme change(end)
    %     calcph_options.ph_mint = .04;       % minimum phoneme duration(sec) bfore looking for a new noise_count
    calcph_options.f_band  = 200;   % +/- Hz around each formant to calculate intensity in sdb for each phoneme
    if use_praat_int_phtrans
        calcph_options.vowel_int_thresh = 66;% - if intensity at midpoint of a phoneme is less than this the phoneme is a consonant
        calcph_options.firstvowel_int_thres = 62; % use this threshold to see if average intensity is low --> consonant
    else
        calcph_options.vowel_int_thresh = 51;% - if intensity at midpoint of a phoneme is less than this the phoneme is a consonant
        calcph_options.firstvowel_int_thres = 48; % use this thre
    end
    
    calcph_options.formant_noise_thresh = 0.5; % if solid fraction in a formant is below this limit then noise=1, else noise = 0
    calcph_options.n_formant = 3; % no. of formants to use
    calcph_options.formant_int_thresh = 6; % if intensity of a formant is above this limit then int=1, else int = 0
    calcph_options.formant_weightsvowl = [2; 3; 2;   1; 2; 1;   1; 2; 1;   1; 1; 1;   1; 1; 1];
    calcph_options.formant_weightscons = [1; 1; 1;   4; 4; 3;   1; 1; 1;   1; 1; 1;   3; 2; 1];
    %                                      mF         sF         eF         int        noise
    calcph_options.ph_match_method = 2; % 1: use angle method to match  2: use weighted difference to match
    calcph_options.pitch_gap_t_thresh = 0.051; % threshold in seconds for gap between start of current phoneme to start of pitch

    [ ph_percent ] = ph_distpercent( ph_vowl, ph_cons, calcph_options );
    [ name_phcode, names ] = Addresbk(  );

    
    %%
    % all_words = {'Call', 'Amol', 'Kumar'};
    % all_words_ph = {[7 29 8], [30 9 34 8], [7 37 9 40]};
    
    fig_list = [2 3; 4 5];
    phseq_try = cell(2,2);
    first_constry = zeros(2,2);
    for word_no = 2:3 % loop over first and last name
        
        int_sdbf_max2 = 8000; % Hz above which ignore sdb data for intensity
        int_sdbf_max1 = 4000; % Hz threshold to compare intensity above and below this to see if there is a plosive or fricative
        
        ints1 = sum(sdb_norm(sdb_f<=int_sdbf_max1,:),1);
        ints2 = sum(sdb_norm((sdb_f<=int_sdbf_max2)&(sdb_f>int_sdbf_max1),:),1);
        
        
        
        %ints2(ints2>=.95*ints1) = 15;%ints2(ints2>=.95*ints1) -10;
        ints = ints1 + ints2;
        %ints = ints1 + min(ints1,ints2);
        figure(7); clf
        plot(sdb_t, [ints1(:) ints2(:) ints(:)]);
        
        %     word_selected = all_words{word_no};
        word_t = formant_t(zstart(word_no):zend(word_no));
        word_ints1 = ints1(zstart(word_no):zend(word_no));
        word_ints2 = ints2(zstart(word_no):zend(word_no));
        thresh = [30 40 50 60]/2.4; % Threshold for noise - higher thresh for higher formants
        word_formants = formants(zstart(word_no):zend(word_no),:);
        word_pitch = pitch(zstart(word_no):zend(word_no),:);
        word_intensity = intensity(zstart(word_no):zend(word_no),:);
        word_intensityslope2 = intensityslope2(zstart(word_no):zend(word_no),:);
        word_intensityslope1 = intensityslope1(zstart(word_no):zend(word_no),:);
        
        [ solid_word, noise_word, mean_word ] = Solid_formant( word_formants, thresh );
        
        figure(fig_list(word_no-1,1)); clf; clear hh1
        noise_plot = noise_word + repmat([0 2 4 6], length(word_t),1);
        for i = 1:ncol
            plot(word_t, word_formants(:,i),linestdot{i})
            hold on;
            plot(word_t, mean_word(:,i),linestsolid{i})
            plot(word_t, solid_word(:,i),linestcircle{i})
        end
        if word_no == 2
            title('Formants for First Name');
        else
            title('Formants for Last Name');
        end
        
        figure(fig_list(word_no-1,2)); clf; clear hh2
        plot(formant_t, [pitch intensity]);
        hold on;
        ax =axis;
        for i=1:length(zstart)
            plot([formant_t(zstart(i)); formant_t(zstart(i))],ax(3:4),'r--','linewidth',1.5);
            plot([formant_t(zend(i)); formant_t(zend(i))],ax(3:4),'r--','linewidth',1.5);
        end
        xlabel('Time (s)'); ylabel('Freq (Hz)'); title('Intensity, Pitch & Word Segments');
        legend('pitch','int','Location','Best')
        if word_no == 2
            title('Pitch & Intensity for First Name');
        else
            title('Pitch & Intensity for Last Name');
        end

        
%         ints = sum(sdb_norm(sdb_f<=int_sdbf_max,:),1);
%         ints  = int_wt_freq_sdbf'*sdb_norm;
        
        if ~use_praat_int_phtrans
            alpha = 0.15; %0.15; %0.15;
            a = [1 -(1-alpha)]; b = alpha; int_filt = filtfilt(b,a,ints);
            int_filt_word = interp1(sdb_t, int_filt, word_t);
        else
            int_filt_word = word_intensity;
        end

        
%             % use filtered intesity for slopes
%             word_intensityslope2 = [0; 0; diff(diff(int_filt_word))]; % used to check first_consonant if pitch gap is too small
%             word_intensityslope1 = [0; diff(int_filt_word)]; % used to check first_consonant if pitch gap is too small
%         
        
        figure(fig_list(word_no-1,2)); hold on;
        plot(word_t,int_filt_word,'g-');
        hold off
        axis('auto')
        peak_dist_time = 0.05; % minimum distance (s) between peaks
        peak_dist = round(peak_dist_time/(word_t(2)-word_t(1))); % min distance in samples
        
        if ~use_praat_int_phtrans
            peak_prominence = 1.0; %3.0; % minimum intensity change around peaks
        else
            peak_prominence = 0.4;
        end
        
        max_wordph = 10; % maximum number of phonemes expected
        
        pitch_start = find(word_pitch > 0, 1, 'first');
        pitch_t_gap = word_t(pitch_start) - word_t(1);
        % above this threshold - consonant, below - consonant or vowel
        peakfnd_start = word_t(1);
        first_consonant = 0;
        if pitch_t_gap > calcph_options.pitch_gap_t_thresh
            first_consonant = 1; % if consonant = 1 then it is a consonant, if it is 0 then it is a vowel or a consonant with pitch
        else
            time_gap = word_t(1) + 0.08; %calcph_options.pitch_gap_t_thresh;
            time_gap_ind = find(word_t >= time_gap, 1, 'first');
            %             if int_filt_word(time_gap_ind) < calcph_options.vowel_int_thresh % word_intensity
            if mean(int_filt_word(1:time_gap_ind)) < calcph_options.firstvowel_int_thres % word_intensity
                first_consonant = 1;
            else
                % check if the average intesity at beginning of word is too low, if so then also a consonant
                int_time_interval = 0.015; % time interval in s over which to calculate average intensity
                int_samp_interval = round(int_time_interval/(word_t(2)-word_t(1)));
                word_intensityslope12 = word_intensityslope2.*word_intensityslope1;
                k1 = find(word_intensityslope12(1:pitch_start) < 0.0,1,'first');
                if isempty(k1)
                    first_consonant = 0;
                else
                    if any((word_intensityslope12(k1:time_gap_ind) > 0) & (word_intensityslope2(k1:time_gap_ind) > 0))
                        [max_slop12,max_slop12loc] = max(word_intensityslope12(1:time_gap_ind));
                        int_samp_start = max_slop12loc;
                        if max_slop12 > 0.039
                            first_consonant = 1;
                        else
                            first_consonant = 0;
                        end
                        
                    else
                        int_samp_start= pitch_start;
                        first_consonant = 0;
                    end
                    %             if mean(word_intensity(int_samp_start+(1:int_samp_interval))) < calcph_options.vowel_int_thresh
                    %                 first_consonant = 1;
                    %             end
                end
            end
        end
        
        
        word_matchpct = zeros(2,1);
        
        
        for firstcons_try = 1:2
            if first_consonant == 1
                calcph_options.prev_ph_type = 1; % trick word_phoneme to think prev was vowel so current is cons
            else
                calcph_options.prev_ph_type = 2; % trick word_phoneme to think prev was cons so current is vowel
            end
            [ ph_trans, ph_start, locs1, locs2, ph_trans_ind ] = Calc_ph_trans( int_filt_word, word_t, peakfnd_start, peak_dist, peak_prominence, first_consonant, pitch_start ); %int_filt_word
            [ ph_trans, force_frstcons,fricplo ] = new_phtrans_FricPlo( ph_trans, ph_trans_ind, word_ints1, word_ints2, first_consonant );
            if( force_frstcons == 1)
                calcph_options.prev_ph_type = 1;
                first_consonant = 1;
            end
            figure(fig_list(word_no-1,1));
            if exist('hh1','var')
                delete(hh1);
            end
            hh1 = add_ylines(ph_trans,'k--');
            hold on; plot(word_t,int_filt_word,'g-'); hold off
            xlabel('Time(s)'); ylabel('Freq (Hz)');% title(['Phoneme segments for the word "' word_selected '"']);
            
            figure(fig_list(word_no-1,2));
            if exist('hh2','var');
                delete(hh2);
            end
            hh2 = add_ylines(ph_trans,'k--');
            
            word_phseq = zeros(max_wordph,1);
            word_phend = zeros(max_wordph,1);
            word_phtype = zeros(max_wordph,1);
            word_phinfo = zeros(max_wordph,20);
            word_phi = 1;
            
            for i = 1:length(ph_trans)
                if i == 1
                    ph_t_start = word_t(1);
                    ph_t_end = ph_trans(1);
                else
                    ph_t_start = ph_trans(i-1);
                    ph_t_end = ph_trans(i);
                end
                
                [ ph_type, phoneme, phoneme_info ] = word_phoneme(word_t, mean_word, solid_word, word_pitch, word_intensity, ph_t_start, ph_t_end, calcph_options, sdb_norm, formants, sdb_f, sdb_t, ph_cons, ph_vowl, ph_vowl_norm, ph_vowl_scale, ph_cons_norm, ph_cons_scale,fricplo(i));
                calcph_options.prev_ph_type = ph_type;
                word_phseq(word_phi) = phoneme;
                word_phend(word_phi) = ph_t_end;
                word_phtype(word_phi) = ph_type;
                word_phinfo(word_phi,1:length(phoneme_info)) = phoneme_info';
                word_phi = word_phi + 1;
            end
            %word_phseq(1:(word_phi - 1))'
            
            phseq_try{firstcons_try, (word_no-1)} = word_phseq(1:(word_phi - 1));
            first_constry(firstcons_try, (word_no-1)) = first_consonant;
%             if word_no == 2
%                 first_flag = 1;
%             else
%                 first_flag = 0;
%             end
%             word_matchpct(firstcons_try) = best_match_addresbk_firstlast( phseq_try{firstcons_try}, name_phcode, ph_percent, first_flag );
            
            first_consonant = ~first_consonant;
        end
        
%         if word_matchpct(2) > (word_matchpct(1)+5)
%             if word_no == 2
%                 found_ph.first = phseq_try{2};
%             elseif word_no == 3
%                 found_ph.last = phseq_try{2};
%             end
%             first_cons(inames,(word_no-1)) = ~first_consonant;
%         else
%             if word_no == 2
%                 found_ph.first = phseq_try{1};
%             elseif word_no == 3
%                 found_ph.last = phseq_try{1};
%             end
%             first_cons(inames,(word_no-1)) = first_consonant;
%         
%         end
        
    end
    
    %%
    prev_bestboth = -1000;
    bestindboth   = 0;
    for firstname_cons_try = 1:2
        for lastname_cons_try = 1:2
            found_ph_try.first = phseq_try{firstname_cons_try,1};
            found_ph_try.last = phseq_try{lastname_cons_try,2};
            [ match_percent_try ] = best_match_addresbk( found_ph_try, name_phcode, ph_percent);
            [best_bothname, bestind] = max(match_percent_try.both);
            if best_bothname > prev_bestboth
               match_percent = match_percent_try;
               found_ph = found_ph_try;
               first_percent = match_percent.first;
               last_percent = match_percent.last;   
               both_percent = match_percent.both;
               prev_bestboth = match_percent.both(bestind);
               bestindboth = bestind;
               best_match_percent(inames,:) = [best_bothname first_percent(bestind) last_percent(bestind)];
               best_match_name{inames} = [names.first{bestind} ' ' names.last{bestind}];
               first_name_firstcons(inames) = first_constry(firstname_cons_try, 1);
               last_name_firstcons(inames) = first_constry(lastname_cons_try, 2);
            end
            
        end
    end
    
    show_match = [(1:length(first_percent))' first_percent last_percent both_percent];
    
    res_str = {};
    res_str{1,1} = sprintf(['\t%-15s','   ','%31s'] ,'  NAME  ', 'First%      Last%       Both%');
    disp(res_str{1,1});
    for iaddr = 1:size(show_match,1)
        res_str{iaddr+1,1} = sprintf(['\t%-15s',' : ','%31s'] , [names.first{iaddr} ' ' names.last{iaddr}], num2str(show_match(iaddr,2:end),5));
        if iaddr == bestindboth
            res_str{iaddr+1,1} = [res_str{iaddr+1,1} ' !!!'];
        end
        disp(res_str{iaddr+1,1});
%         res_str{iaddr,1} = sprintf('%s \t %5.2f \t %5.2f \t %5.2f',[names.first{iaddr} ' ' names.last{iaddr}],show_match(iaddr,2), show_match(iaddr,3),show_match(iaddr,4));
    end    
    
    
    fprintf('Best match percent : %4.1f\n', best_match_percent(inames,1))
    fprintf('Best match name : %s \n', best_match_name{inames})
    disp(['Tested with ' wav_fname]);
    disp('-----------');
    %disp('Correct!');
    
    catch
        error_flag(inames) = 1;
        disp(lasterr);
    end
end

% time_taken = toc

%%
if wav_file_use == 0
    xlswrite('Contact_Names_New.xlsx',correct_name,'Results','A2');
    xlswrite('Contact_Names_New.xlsx',{'Correct Name/file'},'Results','A1');
    xlswrite('Contact_Names_New.xlsx',best_match_name,'Results','B2');
    xlswrite('Contact_Names_New.xlsx',{'BestMatch Name'},'Results','B1');
    xlswrite('Contact_Names_New.xlsx',best_match_percent,'Results','C2');
    xlswrite('Contact_Names_New.xlsx',{'Both%', 'First%', 'Last%'},'Results','C1');
    xlswrite('Contact_Names_New.xlsx',[first_name_firstcons(:) last_name_firstcons(:)],'Results','F2');
    xlswrite('Contact_Names_New.xlsx',{'First cons', 'Last cons'},'Results','F1');
end
