function Ug = gen_fourier(fs,L)
    %clear all

    %fs = 44100;
    %L = 1024;

    first = (1:L/2);
    first_mirror = L/2:-1:1;
    last = (L/2+1:L);

    Ug_fft = zeros(L,1);


    for i = first
        Ug_fft(i) = i;
    end


    Ug_fft(last) = Ug_fft(first_mirror);

    % add phase

%    for i = 1:L
 %      Ug_fft(i) = Ug_fft(i)*exp(2*pi*1i*i/L);
  %  end

    Ug = circshift(ifft(Ug_fft,L,'symmetric'),L/2);

   % plot(Ug)
end