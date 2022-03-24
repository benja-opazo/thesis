% Generates a glottal impulse given an Rd value
function [G, G_time, G_time_int] = GlottalTest(fs, Rd, T0, Ee, f)
    %fs = 48000;
    %Rd = 1;
    %T0 = 1;
    %Ee = 1;
    %f = 1:512;

    % Calculates the T0 normalized te, tp and ta parameters
    [te, tp, ta] = Rd2tetpta(Rd);


    % Calculates the LF model glottal impulse on the frequency domain
    G = gfm_spec_lf(f+1, fs, T0, Ee, te, tp, ta);
    G_int = G./(2*pi*i*f);
    
    % Calculates the glottal impulse on the time domain
    G_shift = [G G(end: -1 : 2)];
    G_int_shift = [G_int G_int(end: -1 : 2)];
    %plot(abs(G_shift))

    G_time = ifft(G_shift,'symmetric')';
    G_time_int_num = cumtrapz(T0/150,G_time);
    G_time_int = ifft(G_int_shift,'symmetric')';
    G = G';
    
    %% Plots
    
    figure(1)
    clf
    plot(abs(G))
    hold on
    plot(abs(G_int))
    title('Spectrum')
    legend('G','G int')
    
    figure(2)
    clf
    plot(G_time)
    hold on
    plot(G_time_int)
    plot(G_time_int_num)
    title('Signal')
    legend('G','G int','G int num')
    
    figure(3)
    plot(G_time_int - G_time_int_num)
    title('Difference between numeric and dft integration')
    
    figure(4)
    clf
    plot(isnan(G_time))
    hold on
    plot(isnan(G_time_int))
    plot(isnan(G_time_int_num))
    title('isnan signal')
    legend('G','G int','G int num')
    
    %{
    figure(1)
    plot(G_time_int)
    title('G time int')
    
    figure(2)
    plot(G_time_int_num)
    title('G time int num')
    
    %}
    
    
    %plot(G_time)

    % The first sample is 0, and the pulse is in a variable that is 1024
    % samples long
    %G_rz = G_time - G_time(1);
    %G_1024 = zeros(1024,1);
    %G_1024((end/2-max(f)):(end/2+max(f))) = [0;G_rz;0];


    %% Falta dividir la señal en las bandas y ver como suena
end