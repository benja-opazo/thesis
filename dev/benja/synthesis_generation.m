%% Add path
addpath('../WORLD')
addpath('../python-PESQ')
addpath(genpath('../covarep'))

%% Variables
audio_names = ["LA9011_ENSS"; "LA9023_ENSS"];
vq_names = ["modal", "breathy", "vocal_fry", "disphonia", "rough"];

fs = 44100;

%% Voice Quality Struct

male_vq = struct();

male_vq.modal = struct(     'rd_param',0.35, ...
                            'spectrum_filtering_exponential', 8.5, 'spectrum_filtering_samples', 45, ...
                            'rpp_multiplier_te', 0.95, 'rpp_multiplier_tp', 0.94, 'rpp_multiplier_ta', 1, 'rpp_k', 0.9);
              
male_vq.breathy = struct(   'rd_param', 2.8, ...
                            'spectrum_filtering_exponential', 8.5, 'spectrum_filtering_samples', 45, ...
                            'rpp_multiplier_te', 0.95, 'rpp_multiplier_tp', 0.94, 'rpp_multiplier_ta', 1, 'rpp_k', 0.5, ...
                            'band_aperiodicity_multiplier', 0.5.*[1 1 1 1 1 1 1]');

male_vq.vocal_fry = struct( 'rd_param', 0.35, ...
                            'f0_multiplier', 0.5, ...
                            'spectrum_filtering_exponential', 8.5, 'spectrum_filtering_samples', 45, ...
                            'rpp_multiplier_te', 0.95, 'rpp_multiplier_tp', 0.94, 'rpp_multiplier_ta', 1, 'rpp_k', 0.9, ...
                            'jitter_amplitude', 12, 'jitter_frequency', 50*0.7, ...
                            'shimmer_amplitude', 12, 'shimmer_frequency', 50*0.7);
                       
male_vq.disphonia = struct( 'rd_param', 1.35, ...
                            'spectrum_filtering_exponential', 8.5,  'spectrum_filtering_samples', 45, ...
                            'rpp_multiplier_te', 0.95,  'rpp_multiplier_tp', 0.94,  'rpp_multiplier_ta', 1,  'rpp_k', 0.5, ...
                            'band_aperiodicity_multiplier', 0.05.*[0.1 1 1 1 1 1 1]');
                       
male_vq.rough = struct(     'rd_param', 0.35, ...
                            'f0_multiplier', 0.95, ...
                            'spectrum_filtering_exponential', 8.5,  'spectrum_filtering_samples', 45, ...
                            'rpp_multiplier_te', 0.95,  'rpp_multiplier_tp', 0.94,  'rpp_multiplier_ta', 1,  'rpp_k', 0.9, ...
                            'jitter_amplitude', 14,  'jitter_frequency', 50*0.5, ...
                            'shimmer_amplitude', 15,  'shimmer_frequency', 50*0.8);
                        
female_vq = struct();

female_vq.modal = struct(   'rd_param', 1, ...
                            'spectrum_filtering_exponential', 12, 'spectrum_filtering_samples', 75, ...
                            'rpp_multiplier_te', 1, 'rpp_multiplier_tp', 1, 'rpp_multiplier_ta', 1, 'rpp_k', 2.5);
                        
female_vq.breathy = struct(     'rd_param', 4, ...
                                'spectrum_filtering_exponential', 12, 'spectrum_filtering_samples', 75, ...
                                'rpp_multiplier_te', 0.8, 'rpp_multiplier_tp', 1, 'rpp_multiplier_ta', 1, 'rpp_k', 0.8, ...
                                'band_aperiodicity_multiplier', 0.5.*[1 1 1 1 1 1 1]');

female_vq.vocal_fry = struct(   'rd_param', 1, ...
                                'spectrum_filtering_exponential', 12, 'spectrum_filtering_samples', 75, ...
                                'f0_multiplier', 0.35, ...
                                'jitter_amplitude', 18, 'jitter_frequency', 50*0.7, ...
                                'shimmer_amplitude', 15, 'shimmer_frequency', 50*0.7);


female_vq.disphonia = struct(   'rd_param', 1.35, ...
                                'spectrum_filtering_exponential', 10, 'spectrum_filtering_samples', 75, ...
                                'rpp_multiplier_te', 0.45, 'rpp_multiplier_tp', 1, 'rpp_multiplier_ta', 0.1, 'rpp_k', 0.5, ...
                                'band_aperiodicity_multiplier', 0.05.*[0.5 1 1 1 1 1 1]');

female_vq.rough = struct(       'rd_param', 0.35, ...
                                'f0_multiplier', 0.92, ...
                                'spectrum_filtering_exponential', 8, 'spectrum_filtering_samples', 75, ...
                                'jitter_amplitude', 17, 'jitter_frequency', 50*0.8, ...
                                'shimmer_amplitude', 17, 'shimmer_frequency', 50*0.8, ...
                                'rpp_multiplier_te', 0.45, 'rpp_multiplier_tp', 1, 'rpp_multiplier_ta', 0, 'rpp_k', 0.8);

% This variable creates a struct that hass the vq struct of all input
% filenames
all_vq = struct(convertStringsToChars(audio_names(1)), male_vq, convertStringsToChars(audio_names(2)), female_vq);

%% WORLD Options
dio_options = struct('f0_floor', 80,...
                     'f0_ceil', 400,...
                     'target_fs', fs,...
                     'frame_period', 2,... % Time in [ms] between 2 filters
                     'channels_in_octave', 24);
                 
cheaptrick_options = struct('q1', -0.15,...
                            'fft_size', 2048);
                        
D4C_options = struct('threshold', 0.85,...
                     'fft_size', 4096);
                 
seeds_options = struct('fft_size', 1024 * 2 ^ ceil(log2(fs / 48000)),...
                     'number_of_samples', 2 ^ ceil(log2(fs / 2)));

%%
% Synthesis

%audiowrite('../outputs/tmp/x.wav',x./max(abs(x)),fs)
%audiowrite('../outputs/tmp/y.wav',y_d./max(abs(y_d)),fs)
tic
for audio = audio_names'
    fprintf('Analyzing: ' + audio + '\n')
    [x, fs] = audioread('../inputs/' + audio + '.wav');
    
    % DIO
    f0_parameter = Dio(x, fs, dio_options);

    % StoneMask is used to improve DIO
    f0_parameter.f0 = StoneMask(x, fs, f0_parameter.temporal_positions, f0_parameter.f0);

    % CheapTrick (Spectrum Estimation)
    spectrum_parameter = CheapTrick(x, fs, f0_parameter, cheaptrick_options);

    % D4C (Aperiodicity Estimation)
    source_parameter = D4CRequiem(x, fs, f0_parameter, D4C_options);

    % Seeds
    seeds_signals = GetSeedsSignals(fs, seeds_options);
    
    fprintf('Voice Quality: ')
    for vq = vq_names
        fprintf(vq)
        
        vq_options = all_vq.(audio).(vq);
        
        y = SynthesisRequiem(source_parameter, spectrum_parameter, seeds_signals, vq_options, 'off');
        
        audiowrite(convertStringsToChars('../outputs/synthesized/experiment2/' + audio + '_' + vq + '.wav'), y./max(abs(y)),fs)
        fprintf(', ')
    end
    fprintf('\n')
end
toc