L = 1024;
K = 1;
T0 = 1/100;
Rd = 2.7;
fs = 44100;

[Te, Tp, Ta] = Rd2tetpta(Rd);
%Ug = rpp(fs, L, T0, K, Te*T0, 1.5*Tp*T0, 0*Ta*T0, 'off'); %Super inestable
%Ug = rpp(fs, L, T0, K, Te*T0, 2*Te*T0/3, 0); % Div por 0 (fixed)
%Ug = rpp(fs, L, T0, K, Te*T0, 2*Te*T0/3 - 0.0001, 0); % Div por 0 (fixed)
%Ug = rpp(fs, L, T0, K, 0.0023, 0.0025, 0*0.05*T0); % Pulso negativo (fixed)
Ug = rpp(fs, L, T0, 1, 0.2350*T0, 0.002, 0.15*T0);

% Zeros
r_z = 1;
alpha_z = 0;

% Poles
r_p = 0.99;
alpha_p = 0;

b = [1, -2*r_z*cos(alpha_z), r_z.^2];
a = [1, -2*r_p*cos(alpha_p), r_p.^2];

%[b,a] = butter(1,100/(fs/2),'high');
    
close all

figure(1)
plot(Ug)
hold on
plot(filter(b,a,Ug))

figure(2)
freqz(b,a)




%bode(b,a)



%plot(Ug)