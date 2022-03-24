% Time domain implementation of the Rosenberg++ Glottal Model
%
% Description
%  This function generates a glottal pulse on the time domain based on the
%  Rosenberg++ Glottal Model [1], using the GFM discussed on [2]
%   
% Input
%  fs  : [Hz] Sampling frequency.
%  L   : [sample] Number of samples of the output signal
%  K   : Amplitude Coefficient.
%  T0  : [s] Fundamental period.
%  Te  : [s] Minimum of the glottal flow derivative waveform (excitation
%  instant). Also known as Glottal Closure Instant
%  Tp  : Maximum of the glottal flow waveform
%  Ta  : [s] Time constant for the return phase
%
% Output
%  Ug   : R++ time-domain Glottal Flow
%
% References
%  [1] R. Veldhuis "A computationally efficient alternative for the
%      Liljencrants–Fantmodel and its perceptual evaluation", The Journal 
%      of the Acoustical Society of America, 103(1), 566-71, 1998.
%  [2] B. Doval, C. d'Alessandro and N. Henrich, "The spectrum of glottal flow
%      models", Acta acustica united with acustica, 92(6), 1026-1046, 2006.
%
% Author
%  Benjamín Opazo <benjamin.opazo.13@sansano.usm.cl>
%
% TODO
%  -
%

%% Inputs

%clear all
%close all

fs = 48000;
L = 1024; % > T0*fs

%These values are tested!!! do not delete
%{
T0 = 1/110;
K = 0.5;%0.3*10^9; % >0 ??
Te = 0.7*T0; % < T0
Tp = 0.65*Te; % 0.5*Te <= tp <= 3/4*Te*(Te + 4*Ta*D)/(Te + 3*Ta*D)
Ta = 0.5*(T0-Te);%0.1*T0; % > 0
%}
%}
%{
T0 = 1/110;
K = 1; % >0 ??
Te = 0.5*T0; % < T0
Tp = 0.65*Te; % 0.5*Te <= tp <= 3/4*Te*(Te + 4*Ta*D)/(Te + 3*Ta*D)
Ta = 0.001*(T0-Te);%0.1*T0; % > 0
%}

%[K, T0, Te, Tp, Ta] = rpp_param(2, 1/100, 0.65, 0.5, 0.5);
Rd = 1.12;
[Te, Tp, Ta] =  Rd2tetpta(Rd);
T0 = 1;
K = 1;
L = fs;
%}


%% Error checking

if (T0*fs > 2*L) % T0 bounds
   error('T0 cannot be higher than L/fs. Pulse is longer than the output variable')
end

% Factor to calculate Tp bounds
if (0.5*Te > Tp) % Tp bounds [1] eq. 13, left hand inequality
    error('Tp is too low. The Glottal pulse will be negative') 
end

if (Ta ~= 0)
    D = 1 - ((T0 - Te)/Ta)/(exp((T0-Te)/Ta) - 1);
    if Tp > 0.75*Te*(Te + 4*Ta*D)/(Te + 3*Ta*D) % Tp bounds [1] eq. 13, right hand inequality
        %warning('Tp is too high. The Glottal pulse will be too skewed') 
    end
end

if (Te >= T0)
   error('Te must be less than T0') 
end

%% Auxiliary Variables calculations

% Time variable
t = (1/fs):(1/fs):T0;

% Samples until it is reached the minimum of the glottal flow derivative
t_gci = 1:floor(Te*fs);
% Samples from the minimum of the GF derivative until the end of the GF
t_gci_end = (floor(Te*fs)+1):floor(T0*fs);

exp_1 = exp(-(T0 - Te)/Ta);

% Position of the pulse on the output variable
%pulse_position_1 = floor(t_gci + L/2 - length(t)/2); %first half
%pulse_position_2 = floor(t_gci_end + L/2 - length(t)/2); %second half

% Position of the pulse on the output variable. The maximum of the
% pulse (Tp*fs) is at position L/2
pulse_position_1 = floor(t_gci + L/2 - floor(Tp*fs)); %first half
pulse_position_2 = floor(t_gci_end + L/2 - floor(Tp*fs)); %second half


%% Parameters calculation
if (Ta ~= 0)
    D = 1 - ((T0 - Te)/Ta)/(exp((T0 - Te)/Ta) - 1);
    Tx = Te*(1 - (0.5*Te^2 - Te*Tp)/(2*Te^2 - 3*Te*Tp + 6*Ta*(Te-Tp)*D));
else
    Tx = Te*(3*Te - 4*Tp)/(2*(2*Te - 3*Tp));
end

% For some reason, the maximum value of the Glottal Flow is not K, but
% instead Tp^3*(2*Tx - Tp)/3. In order to make the maximum of the Glottal
% Flow depend directly on a variable, this auxiliary variable K_real is
% calculated
K_real = K/(Tp^3*(2*Tx - Tp)/3);
%K_real = K;

%We need to calculate the derivative of Ug at t = Te
Ug_d = 4*K_real*Te*(Tp - Te)*(Tx - Te);

%% Glottal Flow calculation

Ug = zeros(1,L);

Ug(pulse_position_1) = K_real.*t(t_gci).^2.*(t(t_gci).^2 - (4/3).*t(t_gci)*(Tp + Tx) + 2*Tp*Tx);

if (Ta ~= 0)
     Ug(pulse_position_2-1) = Ug(pulse_position_1(end)) + Ta*Ug_d*(1 - exp(-(t(t_gci_end) - Te)/Ta) ...
           - ((t(t_gci_end) - Te)/Ta)*exp_1)/(1 - exp_1);
end

%plot(Ug)
%xlim([0,L])
%plot(t,Ug([pulse_position_1 pulse_position_2]))
plot(Ug)
max(Ug)