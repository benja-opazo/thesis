clear all
%%
test_signal = zeros(10000,1);
test_signal(floor(end/2)) = 1;

syms z;
b = [0 1 -1];
a = [1 1.9 0.9025];
f = tf(b,a,1);

filtered_1 = filter(b,a,test_signal);
filtered_2 = filter([1 -1],[2 1],test_signal);

F1mag = abs(fft(filtered_1));
F1phase = angle(fft(filtered_1));

F2mag = abs(fft(filtered_2));
F2phase = angle(fft(filtered_2));


%% Graficos
figure(1)
pzmap(f)

figure(2)
plot(filtered_1)

figure(3)
plot(filtered_2)