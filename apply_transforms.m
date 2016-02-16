function apply_transforms()
% This function applies flattening, truncating, zero-meaning,
% zero-padding, and windowing (Hamming/Hanning)

  global data;
  global plots;
  global analysis;
  global analysed;

  % Are we getting the FFT of the original timeseries or the flattened timeseries?
  if (analysis.apply_bps)
    analysed.transformed = analysed.flattened;
  else
    analysed.transformed = data.timeseries;
  end % if

  % Are we processing just the visible part?
  if (analysis.only_visible)
    xax = plots.timeseries.axis(1:2);
    min_idx = max([floor(xax(1)/analysed.dt)+1, 1]);
    max_idx = min([floor(xax(2)/analysed.dt)+1, analysed.N]); % <-- Check this
    analysed.transformed = analysed.transformed([min_idx:max_idx]);
  end % if

  % Apply zeromean
  if (analysis.zeromean)
    analysed.transformed = analysed.transformed - mean(analysed.transformed);
  end % if

  % Apply zero-pad
  if (analysis.zeropad)
    n = length(analysed.transformed);
    np = analysis.period / analysed.dt; % Number of bins per pulse
    n_extra_zeros = round(ceil(n/np)*np) - n;
    analysed.transformed = [analysed.transformed; zeros(n_extra_zeros,1)];
  end % if

  % Get the length of the timeseries to be transformed
  n  = length(analysed.transformed);

  % Apply Hamming window
  if (analysis.apply_hamming)
    analysed.transformed = hamming(n) .* analysed.transformed;
  end % if

  % Apply Hanning window
  if (analysis.apply_hanning)
    analysed.transformed = hanning(n) .* analysed.transformed;
  end % if

end % function
