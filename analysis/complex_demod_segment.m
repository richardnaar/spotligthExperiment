function [demodAmp, demodPhase] = complex_demod_segment(x, fs, targetFreq, lpCutoffHz, lpOrder, taperFraction)
% x           : 1 x N vector (single-channel data segment)
% fs          : sampling rate (Hz)
% targetFreq  : frequency to demodulate (Hz)
% lpCutoffHz  : low-pass cutoff for envelope (Hz)
% lpOrder     : Butterworth order
% taperFraction : fraction of samples tapered at each edge (0–0.5)

    x = x(:)';  % ensure row vector
    N = numel(x);
    t = (0:N-1) / fs;  % time vector in seconds

    % --- square-cosine taper (similar to your original compdemod) ---
    if taperFraction > 0
        M2 = round(taperFraction * N);
        if 2*M2 >= N
            M2 = floor((N-1)/2);  % guard: keep a flat middle region
        end

        theta = linspace(pi/2, pi, M2);
        squarecos1 = (cos(theta)).^2;

        flatLen = N - 2*M2;
        if flatLen < 0; flatLen = 0; end

        squarecosfunction = [squarecos1, ones(1, flatLen), fliplr(squarecos1)];
    else
        squarecosfunction = ones(1, N);
    end

    % carrier signals
    carrierSin = sin(2*pi*targetFreq .* t);
    carrierCos = cos(2*pi*targetFreq .* t);

    % multiply by taper and carrier
    Xsin = x .* squarecosfunction .* carrierSin;
    Xcos = x .* squarecosfunction .* carrierCos;

    % --- low-pass filter for envelope ---
    % normalized cutoff for butter (0–1; 1 = Nyquist)
    Wn = lpCutoffHz / (fs/2);
    Wn = min(max(Wn, 0.0001), 0.99);  % clamp to legal range

    [B, A] = butter(lpOrder, Wn, 'low');

    XsinF = filtfilt(B, A, Xsin);
    XcosF = filtfilt(B, A, Xcos);

    % amplitude and phase (note factor 2 for analytic representation)
    demodAmp   = 2 * sqrt(XsinF.^2 + XcosF.^2);
    demodPhase = atan2(XsinF, XcosF);   % atan2 is safer than atan(XsinF./XcosF)

end
