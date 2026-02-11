function [COEFFICIENTS, DOPPLER, RANGE, freq, t] = ...
    Range_Doppler_Process(y, sl, ts, PRI, C)

    y  = y(:).';    
    sl = sl(:).';    

    fs = 1/ts;

    PRI_num = round(PRI/ts);
    PRF = 1/PRI;

    sample_num = length(y);
    Trec = sample_num * ts;

    deltaf = 1/Trec;
    pulse_num = round(Trec/PRI);

    t = 0:ts:Trec-ts;
    freq = -fs/2:deltaf:fs/2-deltaf;

    SL_PRI = sl(1:PRI_num);

    yf = fftshift(fft(y));

    if mod(pulse_num, 2)==0
        yff = circshift(yf, pulse_num/2);
        yf_reshape = buffer(yff, pulse_num).';
        SL_PRIf = fftshift(fft(SL_PRI)).';
        SL_PRIff = circshift(SL_PRIf, pulse_num/2);
        BASE = repmat(conj(SL_PRIff), 1, pulse_num);
        COEFFICIENTS = ifft(yf_reshape .* BASE);
    else
        yff = yf;
        yf_reshape = buffer(yff, pulse_num).';
        SL_PRIf = fftshift(fft(SL_PRI)).';
        SL_PRIff = SL_PRIf;
        BASE = repmat(conj(SL_PRIff), 1, pulse_num);
        COEFFICIENTS = ifft(yf_reshape .* BASE);
    end

    RANGE = C/2*(0:ts:PRI-ts);
    DOPPLER = -PRF/2:deltaf:PRF/2-deltaf;
end
