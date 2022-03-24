addpath('../benja')
clear all
%% Input Parameters
fs = 48000;
T0 = 1/100;
L = 1024;

num_im = 0;
%{
for i = 0.25:0.001:2.5
    num_im = num_im + 1;
    Rd = i;
    [Te, Tp, Ta] =  Rd2tetpta(Rd);
    Ug = rpp(fs, L, 1, 1, Te, Tp, Ta);
    g = plot(linspace(0,T0,fs),abs(Ug));
    text(0.7,0.9,strcat('Rd = ',num2str(Rd)))
    saveas(g, strcat('C:\Users\Benja\Desktop\rd\rd_',int2str(num_im),'.jpg'))
end
%}

%Ug = rpp(fs, L, 1, 1, Te, Tp, Ta);
%% Glottal Impulse
Rd = 0.4;
K = 1;
[Te, Tp, Ta] =  Rd2tetpta(Rd);
Ug1 = rpp(fs, L, T0, K, Te*T0, Tp*T0, Ta*T0);
Rd = 1;
[Te, Tp, Ta] =  Rd2tetpta(Rd);
Ug2 = rpp(fs, L, T0, K, Te*T0, Tp*T0, Ta*T0);
Rd = 4.5;
[Te, Tp, Ta] =  Rd2tetpta(Rd);
Ug3 = rpp(fs, L, T0, K, Te*T0, Tp*T0, Ta*T0);

n = linspace(0,L/fs,L);
Ug = [n; Ug1; Ug2; Ug3]';

%plot(Ug)


%%

csvwrite("../outputs/figures/rosenberg_rd04.csv", [n; Ug1]')
csvwrite("../outputs/figures/rosenberg_rd1.csv", [n; Ug2]')
csvwrite("../outputs/figures/rosenberg_rd45.csv", [n; Ug3]')

% Glottal Impulse Spectrum
%G = fft(Ug,L/2 + 1)';