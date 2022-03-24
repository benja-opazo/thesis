% Parameters calculation of the generic GFM for the Rosenberg++ Model
%
% Description
%  This funcion generates the Rosenberg++ specific parameters [1], given the
%  generic parameters of the GFM [2].
%   
% Input
%  E   : The maximum excitation
%  T0  : [s] Fundamental period.
%  Oq  : The open quotient. Range: [0, 1] -> Useful: [0.3, 0.8]
%  am  : The asymmetry coefficient. Range: [0.5, 0.75]
%  Qa  : The return phase quotient. Range: [0, 1]
%
% Output
%  K   : Amplitude Coefficient.
%  T0  : [s] Fundamental period.
%  Te  : [s] Minimum of the glottal flow derivative waveform (excitation
%  instant).
%  Tp  : Maximum of the glottal flow waveform
%  Ta  : [s] Time constant for the return phase
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

function [K, T0, Te, Tp, Ta] = rpp_param(E, T0, Oq, am, Qa)
    %% Error checking
    if (Oq < 0 || Oq > 1) % Oq Bounds
        error('Oq out of bounds! Range is [0, 1]')
    end
    if (Oq < 0.3 || Oq > 0.8)
       warning('It is recommended that Oq is within the range [0.3, 0.8]') 
    end
    if (am < 0.5 || am > 0.75) % am bounds
        error('am out of bounds! Range is [0.5, 0.75]')
    end
    if (Qa < 0 || Qa > 1) % Qa bounds
        error('Qa out of bounds! Range is [0, 1]')
    end
    
    %% Parameters calculation
    Te = Oq*T0;
    Tp = am*Te;
    if (Qa == 0)
        % Parameters calculation of the generic model (Qa = 0 => ta = 0)
        Ta = 0;
        
        % K specific parameters
        Tx = Te*(3*Te - 4*Tp)/(2*(2*Te - 3*Tp));    % See: [2] Appendix 1.2
    else
    % Parameters calculation (Qa ~= 0 => ta ~= 0)
        Ta = Qa*(1-Oq)*T0;
        
        % K specific parameters
        D = 1 - ((T0 - Te)/Ta)/(exp((T0 - Te)/Ta) - 1);
        Tx = Te*(1 - (0.5*Te^2 - Te*Tp)/(2*Te^2 - 3*Te*Tp + 6*Ta*(Te-Tp)*D));
    end
    
    % For some reason, the maximum value of the Glottal Flow is not
    % given by the equations on [2] Appendix 1.2. The following value
    % was obtained by defining E as the real value of the maximum
    % amplitude of the Glottal Flow, and then calculating K
    % accordingly. The equation Tp^3*(2*Tx - Tp)/3 was obtained
    % symbolically
    %K = E/(Tp^3*(2*Tx - Tp)/3);
    K = E;
end