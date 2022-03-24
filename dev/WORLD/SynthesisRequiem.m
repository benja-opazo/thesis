function y = SynthesisRequiem(source_object, filter_object, seeds_signals, vq_options, report_warning)
% Waveform generation from the estimated parameters
% y = SynthesisRequiem(source_object, filter_object, seeds_signals)
%
% Input
%   source_object : F0 and aperiodicity
%   filter_object : spectral envelope
%   seeds_signals : parameters for excitation signal generation
%   vq_options    : Options for voice quality modifications
%
% Output
%   y : synthesized waveform
%
% 2018/04/04: First version

% This piece of code disables warnings if report_warning == 'off'
skip_warning = 0;
if nargin == 5
        if ~strcmp(report_warning,'on') && ~strcmp(report_warning,'off')
           error('Report Warning variable is either "on" or "off"') 
        end
        if strcmp(report_warning,'off')
            skip_warning = 1;
        end
    end

%% Voice Quality Options
fs = filter_object.fs; % Sampling frequency of synthesized signal

rd = 0;
%GF_filter = 
f0_filter_cutoff_frequency = 0;
f0_filter_order  = 1;
f0_multiplier = 1;

jitter_options = struct('fs', fs, ...
                        'JA', 0, ...
                        'JF', 0);
                    
shimmer_options = struct('fs', fs, ...
                         'SA', 0, ...
                         'SF', 0);
                     
vibrato_options = struct('fs', fs, ...
                         'VA', 0, ...
                         'VF', 0);
                     
spectrum_options = struct('exp', 0, ...
                          'n', 0);
                      
rpp_options = struct('te', 1, ...
                     'tp', 1, ...
                     'ta', 1, ...
                     'k', 1);
                 
band_aperiodicity_multiplier = ones(size(source_object.band_aperiodicity,1),1);

if nargin >= 4
    if isfield(vq_options, 'rd_param') == 1
        if (vq_options.rd_param >= 0 && vq_options.rd_param <= 4)
            if (vq_options.rd_param == 0)
                if ~skip_warning
                    warning('Rd parameter is disabled. Using impulse train as excitation signal') 
                end
            elseif (vq_options.rd_param < 0.35)
                error('Rd parameter cannot be lower than 0.35 or higher than 4')
            end
        else
            error('Rd parameter cannot be lower than 0.35 or higher than 4')
        end
      rd = vq_options.rd_param;
    end
    
    if isfield(vq_options, 'f0_filter_frequency') == 1
      if (vq_options.f0_filter_frequency <= 0 || vq_options.f0_filter_frequency > 1)
         error('F0 filter cutoff is defined between 0 and 1') 
      end
      f0_filter_cutoff_frequency = vq_options.f0_filter_frequency;
    end
    if isfield(vq_options, 'f0_filter_order') == 1
      if (vq_options.f0_filter_order <= 0)
         error('F0 filter order is a positive integer') 
      end
      f0_filter_order = vq_options.f0_filter_order;
    end
    
    if isfield(vq_options, 'f0_multiplier') == 1
      if (vq_options.f0_multiplier <= 0)
         error('F0 multiplier can only be greater than 0') 
      end
      f0_multiplier = vq_options.f0_multiplier;
    end
    
    if isfield(vq_options, 'jitter_amplitude') == 1
      if (vq_options.jitter_amplitude < 0)
         error('Jitter amplitude can only be non negative') 
      end
      if (vq_options.jitter_amplitude >= 100)
          if ~skip_warning
            warning('Jitter amplitude is too high!!') 
          end
      end
      
      jitter_options.JA = vq_options.jitter_amplitude;
    end
    
    if isfield(vq_options, 'jitter_frequency') == 1
      if (vq_options.jitter_frequency < 0)
         error('Jitter frequency cutoff can only be non negative') 
      end
      
      jitter_options.JF = vq_options.jitter_frequency;
    end
    
    if isfield(vq_options, 'shimmer_amplitude') == 1
      if (vq_options.shimmer_amplitude < 0)
         error('Shimmer amplitude can only be non negative') 
      end
      if (vq_options.shimmer_amplitude >= 100)
          if ~skip_warning
            warning('Shimmer amplitude is too high!!') 
          end
      end
      
      shimmer_options.SA = vq_options.shimmer_amplitude;
    end
    
    if isfield(vq_options, 'shimmer_frequency') == 1
      if (vq_options.shimmer_frequency < 0)
         error('Shimmer frequency cutoff can only be non negative') 
      end
      
      shimmer_options.SF = vq_options.shimmer_frequency;
    end
    
    if isfield(vq_options, 'vibrato_amplitude') == 1
      if (vq_options.vibrato_amplitude < 0)
         error('Vibrato amplitude can only be non negative') 
      end
      if (vq_options.vibrato_amplitude >= 10)
          if ~skip_warning
            warning('Vibrato amplitude is too high!! (it should be lower than 10)') 
          end
      end
      
      vibrato_options.VA = vq_options.vibrato_amplitude;
    end
   
    if isfield(vq_options, 'vibrato_frequency') == 1
      if (vq_options.vibrato_frequency < 0)
         error('Vibrato frequency cutoff can only be non negative') 
      end
      
      vibrato_options.VF = vq_options.vibrato_frequency;
    end
        
    if isfield(vq_options, 'spectrum_filtering_exponential') == 1
      if (vq_options.spectrum_filtering_exponential < 0)
         error('The exponential value can only be positive') 
      end
      
      spectrum_options.exp = vq_options.spectrum_filtering_exponential;
    end
    
    if isfield(vq_options, 'spectrum_filtering_samples') == 1
      if (vq_options.spectrum_filtering_samples < 0)
         error('The filtered samples have to be more than 0') 
      end
      
      spectrum_options.n = vq_options.spectrum_filtering_samples;
    end
    
    if isfield(vq_options, 'rpp_multiplier_te') == 1
      if (vq_options.rpp_multiplier_te < 0)
         error('The multiplier can only be a positive number') 
      end
      
      rpp_options.te = vq_options.rpp_multiplier_te;
    end
    
    if isfield(vq_options, 'rpp_multiplier_tp') == 1
      if (vq_options.rpp_multiplier_tp < 0)
         error('The multiplier can only be a positive number') 
      end
      
      rpp_options.tp = vq_options.rpp_multiplier_tp;
    end
    
    if isfield(vq_options, 'rpp_multiplier_ta') == 1
      if (vq_options.rpp_multiplier_ta < 0)
         error('The multiplier can only be a positive number') 
      end
      
      rpp_options.ta = vq_options.rpp_multiplier_ta;
    end
    
    if isfield(vq_options, 'rpp_k') == 1
      if (vq_options.rpp_k < 0)
         error('Pulse amplitude can only be a positive number') 
      end
      
      rpp_options.k = vq_options.rpp_k;
    end
    
    if isfield(vq_options, 'band_aperiodicity_multiplier') == 1
      if (size(vq_options.band_aperiodicity_multiplier,1) ~= size(source_object.band_aperiodicity,1) || size(vq_options.band_aperiodicity_multiplier,2) ~= 1)
         error('The band_aperiodicity_multiplier is an array of the form [1 %i]', size(source_object.band_aperiodicity,1)) 
      end
      
      band_aperiodicity_multiplier = vq_options.band_aperiodicity_multiplier;
    end
end

% The following line checks for incompatible options, e.g.:
% jitter_amplitude is present but not jitter_frequency, or
% f0_filter_frequency is present along with jitter options.

% TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% VQ
% Band Aperiodicity Modification
source_object.band_aperiodicity = source_object.band_aperiodicity.*band_aperiodicity_multiplier;

% Spectrum Filtering (Cheaptrick)
if (spectrum_options.exp ~= 0 && spectrum_options.n > 1)
    for i = 1:spectrum_options.n
        filter_object.spectrogram(i,:) =  filter_object.spectrogram(i,:)*exp(i/(spectrum_options.n/spectrum_options.exp) - spectrum_options.exp);
    end
end

% F0 Filtering
if f0_filter_cutoff_frequency > 0 
   [f0_filter_b, f0_filter_a] = butter(f0_filter_order,f0_filter_cutoff_frequency, 'low');
   
   source_object.f0 = filter(f0_filter_b,f0_filter_a, source_object.f0).*source_object.vuv;
end

% F0 Multiplication
source_object.f0 = source_object.f0*f0_multiplier;

% Jitter
if (jitter_options.JA > 0 && jitter_options.JF > 0)
    source_object.f0 = AddJitter(source_object.f0, jitter_options);
else
    if ~skip_warning
        warning('Jitter is disabled')
    end
end

% Vibrato
if (vibrato_options.VA > 0 && vibrato_options.VF > 0)
    vibrato_samples = linspace(0,2*pi,length(source_object.f0));
    vibrato = vibrato_options.VA*sin(vibrato_samples*vibrato_options.VF);
    source_object.f0 = AddVibrato(source_object.f0, source_object.vuv, vibrato);
else
    if ~skip_warning
        warning('Vibrato is disabled')
    end
end



%%

% Generation of the excitation signal
excitation_signal = GetExcitationSignal(source_object.temporal_positions,...
  filter_object.fs, source_object.f0, source_object.vuv,...
  seeds_signals.pulse, seeds_signals.noise, seeds_signals.spectrum, source_object.band_aperiodicity, ...
  rd, shimmer_options, rpp_options, skip_warning);

% Waveform generation based on the overlap-add method
y = GetWaveform(excitation_signal, filter_object.spectrogram,...
  source_object.temporal_positions, source_object.f0, filter_object.fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Thesis functions
function new_f0 = AddJitter(f0, jitter)
    % Generate Brown noise (Bennane2017)
    % Note that the brown noise is generated using a leaky integrator (lowpass)
    % filter to avoid that it wanders off.
    white_noise = wgn(1,length(f0),jitter.JA);
    
    [a,b] = butter(4,jitter.JF/(jitter.fs/2),'low');
    brown_noise = filter(a,b,white_noise);

    new_f0 = f0 + brown_noise;
    % This line fixes the non zero values that should be zero
    new_f0 = (f0~=0).*new_f0;

function y = AddShimmer(x, shimmer)
    % Generate Brown noise (Bennane2017)
    % Note that the brown noise is generated using a leaky integrator (lowpass)
    % filter to avoid that it wanders off.
    white_noise = wgn(length(x),1,shimmer.SA);
    
    [a,b] = butter(4,shimmer.SF/(shimmer.fs/2),'low');
    brown_noise = filter(a,b,white_noise);
    
    % Adds shimmer
    y = x.*(1+ brown_noise);
    
function new_f0 = AddVibrato(f0, vuv, vibrato)
    new_f0 = (vibrato + f0).*vuv;
        
%%
function excitation_signal = GetExcitationSignal(temporal_positions, fs, f0,...
  vuv, pulse_seed, noise_seed, spectrum_seed, band_aperiodicity, ...
  rd, shimmer, rpp_options, skip_warning)
fft_size = size(pulse_seed, 1);
base_index = -fft_size / 2 + 1 : fft_size / 2;
number_of_aperiodicities = size(pulse_seed, 2);

time_axis = temporal_positions(1) : 1 / fs : temporal_positions(end);
periodic_component = zeros(length(time_axis), 1);
aperiodic_component = zeros(length(time_axis), 1);

[pulse_locations_index, interpolated_vuv] = ...
  TimeBaseGeneration(temporal_positions, f0, fs, vuv, time_axis);

% Band-aperiodicity is resampled at sampling frequency of fs Hz
interpolated_aperiodicity = AperiodicityGenration(temporal_positions,...
  band_aperiodicity, time_axis);

% Generation of the aperiodic component
for i = 1 : number_of_aperiodicities
  noise = GenerateNoise(length(aperiodic_component), noise_seed, i);
  aperiodic_component = aperiodic_component +...
    (noise .* interpolated_aperiodicity(i, 1 : length(aperiodic_component))');
end;

% Generation of the periodic component
for i = 1 : length(pulse_locations_index)
  if interpolated_vuv(pulse_locations_index(i)) <= 0.5 ||...
      interpolated_aperiodicity(1, pulse_locations_index(i)) > 0.999
    continue;
  end;
  noise_size = sqrt(max(1, -pulse_locations_index(i) +...
    pulse_locations_index(min(length(pulse_locations_index), i + 1))));
  output_buffer_index =...
    max(1, min(length(time_axis), pulse_locations_index(i) + base_index));

% The expression f0(floor(pulse_locations_index(i)/(length(time_axis)/length(f0)))
% is the corresponding f0 value given a pulse_locations_index(i). To obtain
% this, note that pulse_locations_index is an index of the time_axis, and
% we know that length(f0) is proportional to length(time_axis)
if (floor(pulse_locations_index(i)/(length(time_axis)/length(f0))) > 0)
    current_f0 = f0(floor(pulse_locations_index(i)/(length(time_axis)/length(f0))));    
else
    current_f0 = f0(1);    
end


  response = GetOnePeriodicExcitation(number_of_aperiodicities, pulse_seed, spectrum_seed, ...
    interpolated_aperiodicity(:, pulse_locations_index(i)), noise_size, ...
    current_f0, rd, fs, rpp_options, skip_warning);
  periodic_component(output_buffer_index) =...
    periodic_component(output_buffer_index) + response;
end;

%Shimmer
if (shimmer.SA > 0 && shimmer.SF > 0)
    periodic_component = AddShimmer(periodic_component, shimmer);
else
    if ~skip_warning
        warning('Shimmer is disabled')
    end
end

excitation_signal = periodic_component + aperiodic_component;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function response = GetOnePeriodicExcitation(number_of_aperiodicities,...
  pulse_seed, spectrum_seed, aperiodicity, noise_size, f0, Rd, fs, rpp_options, skip_warning)
if (Rd ~= 0) % Excitation Signal with Rosenberg++ Pulses
    L = 1024;
    % Generate R++ Pulse
    K = rpp_options.k;
    if f0 < 50
        % When unvoiced, defaults to default_f0
        T0 = 1/500;
    else
        T0 = 1/f0;
    end

    [Te, Tp, Ta] = Rd2tetpta(Rd);
    
    if skip_warning
        Ug = rpp(fs, L, T0, K, rpp_options.te*Te*T0, rpp_options.tp*Tp*T0, rpp_options.ta*Ta*T0, 'off');
    else
        Ug = rpp(fs, L, T0, K, rpp_options.te*Te*T0, rpp_options.tp*Tp*T0, rpp_options.ta*Ta*T0, 'on');
    end
    %Ug = rpp(44100, L, T0, K, 0.25*T0, 0.2*T0, 0.05*T0); %vocal fry

    % Calculate Pulse fft
    fft_Ug = fft(Ug, L)';

    % Generate variables placeholders
    fft_filtered_pulse = zeros(L, 7);
    filtered_pulse = zeros(L, 7);

    for i = 1:number_of_aperiodicities
       % FFT of the pulse seeds (this should be calculated previously and added
       % to the input variables of SynthesisRequiem) <- TODO!

       % This is the filtering of the r++ pulse
       fft_filtered_pulse(:,i) = (spectrum_seed(:,i).*fft_Ug)';
       filtered_pulse(:,i) = real(ifft(fft_filtered_pulse(:,i)));
    end
    
    response = zeros(length(filtered_pulse(:, 1)), 1);

    for i = 1 : number_of_aperiodicities
        response = response + filtered_pulse(:, i) * (1 - aperiodicity(i));
    end
else % Excitation signal with impulse train (original signal)
   response = zeros(length(pulse_seed(:, 1)), 1);
    for i = 1 : number_of_aperiodicities
      response = response + pulse_seed(:, i) * (1 - aperiodicity(i));
    end
end

response = response * noise_size; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = GetWaveform(excitation_signal, spectrogram, temporal_positions,...
  f0, fs)
y = zeros(length(excitation_signal), 1);
fft_size = (size(spectrogram, 1) - 1) * 2;
latter_index = fft_size / 2 + 1 : fft_size;

frame_period_sample =...
  round((temporal_positions(2) - temporal_positions(1)) * fs);
win_len = frame_period_sample * 2 - 1;
half_win_len = frame_period_sample - 1;
win = hanning(win_len);

for i = 2 : length(f0) - 2
  origin = (i - 1) * frame_period_sample - half_win_len;
  safe_index = min(length(y), origin : origin + win_len - 1);
    
  tmp = excitation_signal(safe_index) .* win;
  spec = spectrogram(:, i);
  periodic_spectrum = [spec; spec(end - 1 : -1 : 2)];
  
  tmp_cepstrum = real(fft(log(abs(periodic_spectrum)) / 2));
  tmp_complex_cepstrum = zeros(fft_size, 1);
  tmp_complex_cepstrum(latter_index) = tmp_cepstrum(latter_index) * 2;
  tmp_complex_cepstrum(1) = tmp_cepstrum(1);
  
  spectrum = exp(ifft(tmp_complex_cepstrum));
  response = real(ifft(spectrum .* fft(tmp, fft_size)));

  safe_index = min(length(y), origin : origin + fft_size - 1);
  y(safe_index) = y(safe_index) + response;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pulse_locations_index, vuv_interpolated] =...
  TimeBaseGeneration(temporal_positions, f0, fs, vuv, time_axis)
f0_interpolated_raw = ...
  interp1(temporal_positions, f0, time_axis, 'linear', 'extrap');
vuv_interpolated = ...
  interp1(temporal_positions, vuv, time_axis, 'linear', 'extrap');
vuv_interpolated = vuv_interpolated > 0.5;

default_f0 = 500;
f0_interpolated = f0_interpolated_raw .* vuv_interpolated;
f0_interpolated(f0_interpolated == 0) = ...
  f0_interpolated(f0_interpolated == 0) + default_f0;

total_phase = cumsum(2 * pi * f0_interpolated / fs);
%jitter
%F0_length = length(total_phase);
% Sinusoidal
%{
F0_sin_aux = linspace(0,2*pi,F0_length);
%F0_carrier_freq = mean(f0_interpolated)*0.01;
F0_carrier_freq = 20;
F0_carrier = sin(F0_sin_aux*F0_carrier_freq);
total_phase = total_phase + 2*F0_carrier; % add jitter
%}
% Ruidoso
% add white gaussian noise
%total_phase2 = awgn(total_phase, 75, 'measured');

wrap_phase = rem(total_phase, 2 * pi);
pulse_locations = time_axis(abs(diff(wrap_phase)) > pi);
pulse_locations_index = round(pulse_locations * fs) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function multi_aperiodicity = ...
  AperiodicityGenration(temporal_positions, band_aperiodicity, time_axis)
number_of_aperiodicities = size(band_aperiodicity, 1);
multi_aperiodicity = zeros(number_of_aperiodicities, length(time_axis));

for i = 1 : number_of_aperiodicities
  multi_aperiodicity(i, :) = interp1(temporal_positions,...
    10 .^ (band_aperiodicity(i, :) / 10), time_axis, 'linear');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = GenerateNoise(N, noise_seed, frequency_band)
persistent current_index;

if isempty(current_index)
  current_index = zeros(1, size(noise_seed, 2));
end
noise_length = size(noise_seed, 1);

index =...
  rem(current_index(frequency_band) : current_index(frequency_band) + N - 1,...
  noise_length);
n = noise_seed(index + 1, frequency_band);
current_index(frequency_band) = index(end);
