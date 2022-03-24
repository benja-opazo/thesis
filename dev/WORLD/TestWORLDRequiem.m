%%  Test script for WORLD analysis/synthesis with new waveform generator
% 2018/04/04: First version

%% Add path
addpath('../benja')
addpath('../python-PESQ')
addpath(genpath('../covarep'))
%% Audio Inputs

% Peter
%[x, fs] = audioread('../inputs/fem_normal_peter.wav');
%[x, fs] = audioread('../inputs/fem_breathy_peter.wav');
%[x, fs] = audioread('../inputs/fem_rough_peter.wav');

% aaaaaaaaaaaaaaaaaaaa
%[x, fs] = audioread('../inputs/fem_normal_aa.wav');
%[x, fs] = audioread('../inputs/fem_breathy_aa.wav');
%[x, fs] = audioread('../inputs/fem_rough_aa.wav');

% Papa
%[x, fs] = audioread('../inputs/papa_short_mono.wav');

%[x, fs] = audioread('../inputs/male_runn.wav');
[x, fs] = audioread('../inputs/male_sust.wav');
%[x, fs] = audioread('../inputs/female_runn.wav');
%[x, fs] = audioread('../inputs/female_sust.wav');

%[x, fs] = audioread('../inputs/LA9011_ENSS.wav');


%% DIO (F0 Estimation)
% Options
dio_options = struct('f0_floor', 80,...
                     'f0_ceil', 400,...
                     'target_fs', fs,...
                     'frame_period', 2,... % Time in [ms] between 2 filters
                     'channels_in_octave', 24);
% DIO
f0_parameter_d = Dio(x, fs, dio_options);

% StoneMask is used to improve DIO
f0_parameter_d.f0 = StoneMask(x, fs, f0_parameter_d.temporal_positions, f0_parameter_d.f0);

%% CheapTrick (Spectrum Estimation)
% Options
cheaptrick_options = struct('q1', -0.15,...
                            'fft_size', 2048);
% CheapTrick
spectrum_parameter_d = CheapTrick(x, fs, f0_parameter_d, cheaptrick_options);
%% D4C (Aperiodicity Estimation)
% Options
D4C_options = struct('threshold', 0.85,...
                     'fft_size', 4096);


source_parameter_d = D4CRequiem(x, fs, f0_parameter_d, D4C_options);
%% Synthesis
% Options
seeds_options = struct('fft_size', 1024 * 2 ^ ceil(log2(fs / 48000)),...
                     'number_of_samples', 2 ^ ceil(log2(fs / 2)));

% Seeds
seeds_signals_d = GetSeedsSignals(fs, seeds_options);

% Voice Quality modifications
%{
%vq_options = struct();
%vq_options.rd_param = 0.35;
%vq_options.GF_filter = 0.05; %glottal flow filter frequency cuttoff <- FALTA ESTO???
%vq_options.f0_filter_frequency = 0.1;
%vq_options.f0_filter_order = 8; 
%vq_options.f0_multiplier = 0.5;
%vq_options.jitter_amplitude = 18; vq_options.jitter_frequency = fs/16; 
%vq_options.shimmer_amplitude = 11; vq_options.shimmer_frequency = 50; 
%vq_options.jitter_amplitude = 9; vq_options.jitter_frequency = fs/16; % valores fem_normal
%vq_options.shimmer_amplitude = 12; vq_options.shimmer_frequency = 15;   % valores fem_normal
%vq_options.vibrato_amplitude = 5; vq_options.vibrato_frequency = 25;
%vq_options.spectrum_filtering_exponential = 6; vq_options.spectrum_filtering_samples = 55;
%vq_options.rpp_multiplier_te = 0.45 ; vq_options.rpp_multiplier_tp = 1; vq_options.rpp_multiplier_ta = 0.1; vq_options.rpp_k = 1;
%vq_options.band_aperiodicity_multiplier = [1 1 1 1 1 1 1];
%}

%{
%Original WORLD Implementation
vq_options = struct();
%}

%{
% Modal Voice Papa
vq_options = struct();
vq_options.rd_param = 0.35;
vq_options.spectrum_filtering_exponential = 8; vq_options.spectrum_filtering_samples = 45;
vq_options.rpp_multiplier_te = 1; vq_options.rpp_multiplier_tp = 0.5; vq_options.rpp_multiplier_ta = 1; vq_options.rpp_k = 1;
%}

%{
% Modal Male Sust
vq_options = struct();
%vq_options.rd_param = 0.35;
vq_options.rd_param = 4;
vq_options.spectrum_filtering_exponential = 8.5; vq_options.spectrum_filtering_samples = 45;
vq_options.rpp_multiplier_te = 0.95; vq_options.rpp_multiplier_tp = 0.94; vq_options.rpp_multiplier_ta = 1; vq_options.rpp_k = 0.9;
%}

%{
% Breathy Voice Papa
vq_options = struct();
vq_options.rd_param = 2.5;
vq_options.spectrum_filtering_exponential = 10; vq_options.spectrum_filtering_samples = 55;
vq_options.rpp_multiplier_te = 1 ; vq_options.rpp_multiplier_tp = 1; vq_options.rpp_multiplier_ta = 1; vq_options.rpp_k = 1;
%}

%
% Breathy Male Sust
vq_options = struct();
vq_options.rd_param = 2.8;
vq_options.spectrum_filtering_exponential = 8.5; vq_options.spectrum_filtering_samples = 45;
vq_options.rpp_multiplier_te = 0.95; vq_options.rpp_multiplier_tp = 0.94; vq_options.rpp_multiplier_ta = 1; vq_options.rpp_k = 0.5;
vq_options.band_aperiodicity_multiplier = 0.5.*[1 1 1 1 1 1 1]';
%}


%{
% Vocal Fry Papa
vq_options = struct();
vq_options.rd_param = 0.35;
vq_options.f0_multiplier = 0.5;
vq_options.jitter_amplitude = 18; vq_options.jitter_frequency = fs/16; 
vq_options.shimmer_amplitude = 11; vq_options.shimmer_frequency = 50; 
vq_options.spectrum_filtering_exponential = 6; vq_options.spectrum_filtering_samples = 55;
%}

%{
% Vocal Fry Male Sust
vq_options = struct();
vq_options.rd_param = 0.35;
vq_options.f0_multiplier = 0.5;
vq_options.spectrum_filtering_exponential = 8.5; vq_options.spectrum_filtering_samples = 45;
vq_options.rpp_multiplier_te = 0.95; vq_options.rpp_multiplier_tp = 0.94; vq_options.rpp_multiplier_ta = 1; vq_options.rpp_k = 0.9;
vq_options.jitter_amplitude = 12; vq_options.jitter_frequency = 50*0.7;
vq_options.shimmer_amplitude = 12; vq_options.shimmer_frequency = 50*0.7;
%}

%{
% Disphonia Papa
vq_options = struct();
vq_options.rd_param = 1.35;
vq_options.spectrum_filtering_exponential = 6; vq_options.spectrum_filtering_samples = 55;
vq_options.rpp_multiplier_te = 0.45 ; vq_options.rpp_multiplier_tp = 1; vq_options.rpp_multiplier_ta = 0.1; vq_options.rpp_k = 0.05;
%}

%{
% Disphonia Male Sust
vq_options = struct();
vq_options.rd_param = 1.35;
vq_options.spectrum_filtering_exponential = 8.5; vq_options.spectrum_filtering_samples = 45;
vq_options.rpp_multiplier_te = 0.95; vq_options.rpp_multiplier_tp = 0.94; vq_options.rpp_multiplier_ta = 1; vq_options.rpp_k = 0.5;
vq_options.band_aperiodicity_multiplier = 0.05.*[0.1 1 1 1 1 1 1]';
%}

%{
% Rough Voice Papa
vq_options = struct();
vq_options.rd_param = 0.35;
vq_options.f0_multiplier = 0.95;
vq_options.spectrum_filtering_exponential = 6; vq_options.spectrum_filtering_samples = 55;
vq_options.jitter_amplitude = 13; vq_options.jitter_frequency = fs/8; 
vq_options.shimmer_amplitude = 13; vq_options.shimmer_frequency = 50; 
vq_options.rpp_multiplier_te = 0.45 ; vq_options.rpp_multiplier_tp = 1; vq_options.rpp_multiplier_ta = 0; vq_options.rpp_k = 0.8;
%}

%{
% Rough Voice Male Sust
vq_options = struct();
vq_options.rd_param = 0.35;
vq_options.f0_multiplier = 0.95;
vq_options.spectrum_filtering_exponential = 8.5; vq_options.spectrum_filtering_samples = 45;
vq_options.rpp_multiplier_te = 0.95; vq_options.rpp_multiplier_tp = 0.94; vq_options.rpp_multiplier_ta = 1; vq_options.rpp_k = 0.9;
vq_options.jitter_amplitude = 14; vq_options.jitter_frequency = 50*0.5;
vq_options.shimmer_amplitude = 15; vq_options.shimmer_frequency = 50*0.8;
%}
% ----------------------------------------------------------------------------------------------------------------------------------
%{
%Modal Voice Fem Normal
vq_options = struct();
vq_options.rd_param = 1;
vq_options.spectrum_filtering_exponential = 12; vq_options.spectrum_filtering_samples = 75;
vq_options.rpp_multiplier_te = 1 ; vq_options.rpp_multiplier_tp = 1; vq_options.rpp_multiplier_ta = 1; vq_options.rpp_k = 2.5;
%}

%{
% Breathy Voice Fem Normal
vq_options = struct();
vq_options.rd_param = 4;
vq_options.spectrum_filtering_exponential = 12; vq_options.spectrum_filtering_samples = 75;
vq_options.rpp_multiplier_te = 0.8 ; vq_options.rpp_multiplier_tp = 1; vq_options.rpp_multiplier_ta = 1; vq_options.rpp_k = 0.8;
vq_options.band_aperiodicity_multiplier = 0.5.*[1 1 1 1 1 1 1]';
%}

%{
% Vocal Fry Fem Normal
vq_options = struct();
vq_options.rd_param = 1;
vq_options.spectrum_filtering_exponential = 12; vq_options.spectrum_filtering_samples = 75;
vq_options.f0_multiplier = 0.35;
vq_options.jitter_amplitude = 18; vq_options.jitter_frequency = 50*0.7; 
vq_options.shimmer_amplitude = 15; vq_options.shimmer_frequency = 50*0.7; 
%}

%{
%Disphonia Fem Normal
vq_options = struct();
vq_options.rd_param = 1.35;
vq_options.spectrum_filtering_exponential = 10; vq_options.spectrum_filtering_samples = 75;
vq_options.rpp_multiplier_te = 0.45 ; vq_options.rpp_multiplier_tp = 1; vq_options.rpp_multiplier_ta = 0.1; vq_options.rpp_k = 0.5;
vq_options.band_aperiodicity_multiplier = 0.05.*[0.5 1 1 1 1 1 1]';
%}

%{
% Rough Voice Fem Normal
vq_options = struct();
vq_options.rd_param = 0.35;
vq_options.f0_multiplier = 0.92;
vq_options.spectrum_filtering_exponential = 8; vq_options.spectrum_filtering_samples = 75;
vq_options.jitter_amplitude = 17; vq_options.jitter_frequency = 50*0.8; 
vq_options.shimmer_amplitude = 17; vq_options.shimmer_frequency = 50*0.8;
vq_options.rpp_multiplier_te = 0.45 ; vq_options.rpp_multiplier_tp = 1; vq_options.rpp_multiplier_ta = 0; vq_options.rpp_k = 0.8;
%}

% Synthesis
y_d = SynthesisRequiem(source_parameter_d, spectrum_parameter_d, seeds_signals_d, vq_options, 'off');

%soundsc(y_d,fs)
%soundsc([x./max(abs(x)); y_d./max(abs(y_d))],fs)

%audiowrite('../outputs/tmp/x.wav',x./max(abs(x)),fs)
%audiowrite('../outputs/tmp/y.wav',y_d./max(abs(y_d)),fs)

%filename = char(strcat('../outputs/synthesized/', string(datetime('now', 'Format', 'yyyyMMdd-hhmmss')), '.wav'));
%audiowrite(filename, y_d./max(abs(y_d)), fs)

%% Spectrum parameter
%{
plot(spectrum_parameter_d.spectrogram(:,100))
hold on
plot(spectrum_parameter_d2.spectrogram(:,100))
%}
%% CPP
%{
% CPP with smoothing
cpp_x = cpp(x,fs,1);
cpp_yd = cpp(y_d,fs,1);

figure
hold on
plot(cpp_x(:,2),cpp_x(:,1) - 0*mean(cpp_x(:,1)))
plot(cpp_yd(:,2),cpp_yd(:,1) - 0*mean(cpp_yd(:,1)))
legend('Input','Output')
title('CPP of voice segment')
xlabel('Time [samples]')
ylabel('CPP')
%}

%% FFT Comparison
%{
figure
plot(10*log10(abs(fft(x))/max(abs(fft(x)))))
xlim([0, length(x)/2])
title('FFT input')

figure
plot(10*log10(abs(fft(y_d))/max(abs(fft(y_d)))))
xlim([0, length(y_d)/2])
title('FFT synthesized')
%}

%%
%
x_fft = 10*log10(abs(fft(x))/max(abs(fft(x))));
y_d_fft = 10*log10(abs(fft(y_d))/max(abs(fft(y_d))));

x_fft = x_fft(1:end/2);
y_d_fft = y_d_fft(1:end/2);

x_f = linspace(0,fs/2,length(x_fft))';
y_d_f = linspace(0,fs/2,length(y_d_fft))';

%csvwrite("../outputs/figures/x_fft.csv", [x_f x_fft])
csvwrite("../outputs/figures/y_breathy_fft.csv", [y_d_f y_d_fft])
%}
%% Waveforms Comparison
%{
figure
plot(x)
title('Waveform of x')

figure
plot(y_d)
title('Waveform of yd')
%}

%% Rd parameter test
%{
%https://tel.archives-ouvertes.fr/tel-00554763/document (fig 2.5)
[Te, Tp, Ta] = Rd2tetpta(0.3);
T0 = 1/100;
Ug = rpp(48000, 1024, T0, 100, Te*T0, Tp*T0, Ta*T0);

figure
plot(diff(Ug))
title("Rd = 0.3")

[Te, Tp, Ta] = Rd2tetpta(1);
Ug = rpp(48000, 1024, 1/100, 100, Te*T0, Tp*T0, Ta*T0);

figure
plot(diff(Ug))
title("Rd = 1")

[Te, Tp, Ta] = Rd2tetpta(2.7);
Ug = rpp(48000, 1024, 1/100, 100, Te*T0, Tp*T0, Ta*T0);

figure
plot(diff(Ug))
title("Rd = 2.7")
%}