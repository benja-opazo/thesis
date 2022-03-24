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
% 2020/11/17: First version
%
% TODO
%  -
%

%% Inputs
%{
clear all
%close all

fs = 48000;
L = 1024;
T0 = 1/100;
t = (1/fs):(1/fs):T0;
K = 1;
Te = T0/2; % < T0
Tp = T0/2;
Ta = 0;
%}
function Ug = rpp(fs, L, T0, K, Te, Tp, Ta, report_warning)
    %% This piece of code disables warnings if report_warning == 'off'
    skip_warning = 0;
    if nargin == 8
        if ~strcmp(report_warning,'on') && ~strcmp(report_warning,'off')
           error('Report Warning variable is either "on" or "off"') 
        end
        if strcmp(report_warning,'off')
            skip_warning = 1;
        end
    end

    %% Error checking
    if (T0*fs > 2*L) % T0 bounds
       error('T0 cannot be higher than 2*L/fs. Pulse is longer than the output variable')
    end

    % Factor to calculate Tp bounds
    if (0.5*Te > Tp) % Tp bounds [1] eq. 13, left hand inequality
        error('Tp is too low. The Glottal pulse will be negative. To fix this, Tp should be higher than %f (Te/2)', Te/2) 
    end

    if (Ta ~= 0)
        D = 1 - ((T0 - Te)/Ta)/(exp((T0-Te)/Ta) - 1);
        %D = 1 - ((T0 - Te)/Ta)/expm((T0-Te)/Ta);
        if Tp > 0.75*Te*(Te + 4*Ta*D)/(Te + 3*Ta*D) % Tp bounds [1] eq. 13, right hand inequality
            if ~skip_warning
                warning('Tp is too high. The Glottal pulse will be too skewed') 
            end
        end
    else
        if (Tp == 2*Te/3) 
           error('If Ta == 0 and Tp == 2Te/3, then the R++ Model transforms into R+ model (see Velhuis 1998, eq. 11). This case generates a division by 0')
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
    % calculated. Note that divisions by 0 can occurr.
    
    if (Ta ~= 0)
        if abs(Tp - 2*Tx) <= 0.001
            if ~skip_warning
                warning('Instability calculating the amplitude of the pulse: division by 0. To solve this, Tp should not be close to %.12f given these T0, Te and Ta parameters', 2*Tx)
            end
        end
        if (Tp < (4*D*Ta*Te + Te^2)/(2*D*Ta+Te))
            K_real = K/(Tp^3*(2*Tx - Tp)/3);
        else
            K_real = -K/(Tp^3*(2*Tx - Tp)/3);
        end
    else
        if abs(Tp - Te) <= 0.001
            if ~skip_warning
                warning('Instability calculating the amplitude of the pulse: division by 0. To solve this, Tp should not be close to %.12f (Te) given these T0, Te and Ta parameters', Te)
            end
        end
        %if (Tp < 2*Te/3)% && Tp < Te)
        if (Tp < Te)
            K_real = K/(Tp^3*(2*Tx - Tp)/3);
        else
            K_real = -K/(Tp^3*(2*Tx - Tp)/3);
        end
    end
    
    %We need to calculate the derivative of Ug at t = Te
    Ug_d = 4*K_real*Te*(Tp - Te)*(Tx - Te);

    %% Glottal Flow calculation

    Ug = zeros(1,L);

    Ug(pulse_position_1) = K_real.*t(t_gci).^2.*(t(t_gci).^2 - (4/3).*t(t_gci)*(Tp + Tx) + 2*Tp*Tx);

    if (Ta ~= 0)
         Ug(pulse_position_2-1) = Ug(pulse_position_1(end)) + Ta*Ug_d*(1 - exp(-(t(t_gci_end) - Te)/Ta) ...
               - ((t(t_gci_end) - Te)/Ta)*exp_1)/(1 - exp_1);
    end
    Ug(pulse_position_2(end):end) = Ug(pulse_position_2(end-1));
end