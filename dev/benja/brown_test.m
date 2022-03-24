close all
clear all

white_noise_short = wgn(1,100,10);
white_noise_long = wgn(1,48000,10);

bad_brown_short = cumsum(white_noise_short);
bad_brown_long = cumsum(white_noise_long);

    
[a,b] = butter(2,0.1,'low');
good_brown_short = filter(a,b,white_noise_short);
good_brown_long = filter(a,b,white_noise_long);



figure
plot(bad_brown_short)
title('Naive Brownian Noise')
set(gcf,'Position',[0 0 1000 562])

figure
plot(bad_brown_long)
title('Naive Brownian Noise')
set(gcf,'Position',[0 0 1000 562])

figure
plot(good_brown_short)
title('Brownian Noise')
set(gcf,'Position',[0 0 1000 562])

figure
plot(good_brown_long)
title('Brownian Noise')
set(gcf,'Position',[0 0 1000 562])
