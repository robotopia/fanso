function calc_fft()

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  % Are we getting the FFT of the original timeseries or the flattened timeseries?
  if (analysis.apply_bps)
    to_be_ffted = analysed.flattened;
  else
    to_be_ffted = data.timeseries;
  end % if

  % Are we processing just the visible part?
  if (analysis.only_visible)
    xax = plots.timeseries.axis(1:2);
    min_idx = max([floor(xax(1)/analysed.dt)+1, 1]);
    max_idx = min([floor(xax(2)/analysed.dt)+1, analysed.N]); % <-- Check this
    to_be_ffted = to_be_ffted([min_idx:max_idx]);
  end % if

  % Apply zeromean
  if (analysis.zeromean)
    to_be_ffted = to_be_ffted - mean(to_be_ffted);
  end % if

  % Apply zero-pad
  if (analysis.zeropad)
    n = length(to_be_ffted);
    np = analysis.period / analysed.dt; % Number of bins per pulse
    n_extra_zeros = round(ceil(n/np)*np) - n;
    to_be_ffted = [to_be_ffted; zeros(n_extra_zeros,1)];
  end % if

  % Get the length of the timeseries to be FFT'd
  n  = length(to_be_ffted);

  % Apply Hamming window
  if (analysis.apply_hamming)
    to_be_ffted = hamming(n) .* to_be_ffted;
  end % if

  % Apply Hanning window
  if (analysis.apply_hanning)
    to_be_ffted = hanning(n) .* to_be_ffted;
  end % if

  % Calculate the values...
  % ...for the FFT abscissa...
  df                       = 1/(analysed.dt*n);
  analysed.spectrum_freqs  = [0:(n-1)] * df;
  % ...and the FFT ordinate
  analysed.spectrum_vals   = fft(to_be_ffted);

end % function

