function SE = SpectralEnvelopeH1H2(x)
    fft_amplitude_db = 10 * log10( abs( fft( x ) ).^2 );
    
    
    
    [max_h1_db, pos_h1] = max(fft_amplitude_db(2:end-1));
    if (pos_h1 > 15)
        max_h2_db = max(fft_amplitude_db(2*pos_h1 - 30:2*pos_h1 + 30));
    else
        max_h2_db = max(fft_amplitude_db(pos_h1 + 5:2*pos_h1 + 30));
    end
    
    SE = max_h1_db - max_h2_db;
end