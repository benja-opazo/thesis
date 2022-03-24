%Rd validation
addpath(genpath('../covarep'))
%% Common variables
fs = 48000;
rd = 1;
T0 = 1/50;
%L = T0*fs/2 + 220;
L = 1024;
[te, tp, ta] = Rd2tetpta(rd);

%% LF
ncyc = 1;
period = (1/T0);
t = 0:1/period:ncyc;

%ug_lf = glotlf(0,t, [T0*te, T0*tp, T0*ta]);
ug_lf = glotlf(0,t, [te*T0, tp*T0, ta*T0]);
[a,b] = butter(2,0.045, 'low');
ug_lf_filtered = filter(a,b,ug_lf);

plot(diff(ug_lf_filtered))

%% Rpp
%ug_rpp = rpp(fs, L, T0, 1, T0*te, T0*tp, T0*ta);

%plot(diff(ug_rpp))
   
%clf

%hold on
%plot(ug_rpp)

%G = glotlf(

%G_fft = [G(end:-1:1) G];

%G_fft2 = [G G(end:-1:1)];

%plot(ifft(G_fft))
%plot(abs(G_fft))


%plot(circshift(abs(ifft(G_fft2)),1024))

%m2 = imresize(ug_rpp, [1 length(t)], 'nearest');
%plot(m2)

