function SE = SpectralEnvelopeH1H2_covarep(x,fs)
    %fft_amplitude_db = 10 * log10( abs( fft( x ) ).^2 );
    % Settings
    F0min = 80; % Minimum F0 set to 80 Hz
    F0max = 500; % Maximum F0 set to 80 Hz
    frame_shift = 10; % Frame shift in ms
    
    [srh_f0,srh_vuv,srh_vuvc,srh_time] = pitch_srh(x,fs,F0min,F0max,frame_shift);
    
    try
        [creak_pp,creak_bin] = detect_creaky_voice(x,fs); % Detect creaky voice
        creak=interp1(creak_bin(:,2),creak_bin(:,1),1:length(x));
        creak(creak<0.5)=0; creak(creak>=0.5)=1;
        do_creak=1;
    catch
        disp('Version or toolboxes do not support neural network object used in creaky voice detection. Creaky detection skipped.')
        creak=zeros(length(x),1);
        creak_pp=zeros(length(x),2);
        creak_pp(:,2)=1:length(x);
        do_creak=0;
    end
    
    se_gci = se_vq(x,fs,median(srh_f0),creak);           % SE-VQ
    [gf_iaif,gfd_iaif] = iaif_ola(x,fs);    % Glottal flow (and derivative) by the IAIF method
    [NAQ,QOQ,H1H2,HRF,PSP] = get_vq_params(gf_iaif,gfd_iaif,fs,se_gci); % Estimate conventional glottal parameters
    SE = mean(H1H2(:,2));
end