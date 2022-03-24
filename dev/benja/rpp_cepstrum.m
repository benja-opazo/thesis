%% Add path
addpath('../benja')
%% Cepstrum Calculation
L = 1024;
fft_size = 1024;
K = 1;
T0 = 1/100;
Rd = 1;
fs = 44100;

[Te, Tp, Ta] = Rd2tetpta(Rd);
Ug = rpp(fs, L, T0, K, Te*T0, Tp*T0, Ta*T0, 'off');

fft_Ug = fftshift(Ug.^2);   % Fourier Transform
ps_Ug = fft_Ug; % Power Spectrum

logps_Ug = log10(ps_Ug);

cepstrum_Ug = ifftshift(logps_Ug);

ps = GetPowerSpectrum(Ug, fs, fft_size, 1/T0);
ps_smoothed = LinearSmoothing(ps, 1/T0, fs, fft_size);

%plot(cepstrum_Ug)
plot(ps./max(abs(ps)))
hold on
plot(ps_Ug./max(abs(ps_Ug)))

function power_spectrum = GetPowerSpectrum(waveform, fs, fft_size, f0)
    power_spectrum = abs(fft(waveform(:), fft_size)) .^ 2;
    %power_spectrum = power_spectrum;
    % DC correction
    frequency_axis = (0 : fft_size - 1)' / fft_size * fs;
    low_frequency_axis = frequency_axis(frequency_axis <  f0 + fs / fft_size);
    low_frequency_replica = interp1(f0 - low_frequency_axis,...
      power_spectrum(frequency_axis < f0 + fs / fft_size),...
      low_frequency_axis(:), 'linear', 'extrap');
    power_spectrum(frequency_axis < f0) =...
      low_frequency_replica(frequency_axis < f0) +...
      power_spectrum(frequency_axis < f0);
    power_spectrum(end : -1 : fft_size / 2 + 2) = power_spectrum(2 : fft_size / 2);
end


function smoothed_spectrum = LinearSmoothing(power_spectrum, f0, fs, fft_size)
    double_frequency_axis = (0 : 2 * fft_size - 1)' / fft_size * fs - fs;
    double_spectrum = [power_spectrum; power_spectrum];

    double_segment = cumsum(double_spectrum * (fs / fft_size));
    center_frequency = (0 : fft_size / 2)' / fft_size * fs;
    low_levels = interp1H(double_frequency_axis + fs / fft_size / 2,...
      double_segment, center_frequency - f0 / 3);
    high_levels = interp1H(double_frequency_axis + fs / fft_size / 2,...
      double_segment, center_frequency + f0 / 3);

    smoothed_spectrum = (high_levels - low_levels) * 1.5 / f0;
    smoothed_spectrum =...
      smoothed_spectrum + abs(randn(length(smoothed_spectrum), 1)) * eps;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the imprementation of a matlab function
function yi = interp1H(x, y, xi)
    delta_x = x(2) - x(1);
    xi = max(x(1), min(x(end), xi));
    xi_base = floor((xi - x(1)) / delta_x);
    xi_fraction = (xi - x(1)) / delta_x - xi_base;
    delta_y = [diff(y); 0];
    yi = y(xi_base + 1) + delta_y(xi_base + 1) .* xi_fraction;
end