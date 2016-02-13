function flatten()

  global data;
  global analysis;
  global analysed;

  if (isempty(analysis.breakpoints))
    analysed.flattened = data.timeseries;
    analysed.ms = [];
    analysed.cs = [];
    return
  end % if

  % Set up flattened vector
  analysed.flattened = zeros(size(data.timeseries));

  % Book-end the breakpoints with initial and final values
  t   = analysed.t;
  bps = [t(1)-1, analysis.breakpoints, t(end)+1];

  % Calculate the breakpoint mask from the profile mask
  ph = mod(t, analysis.period) / analysis.period;
  pm = analysis.profile_mask;
  if (~isempty(pm))
    if (pm(1) <= pm(2))
      mask_idxs = (ph >= pm(1)) & (ph <= pm(2));
    else
      mask_idxs = (ph >= pm(1)) | (ph <= pm(2));
    end % if
    breakpoint_mask = ~mask_idxs;
  end % if

  % Loop through the breakpoints and detrend each segment
  ms = zeros(length(bps)-1, 1);
  cs = zeros(length(bps)-1, 1);
  for i = 1:(length(bps)-1)
    flatten_idxs  = (t >= bps(i)) & (t < bps(i+1));
    model_idxs    = flatten_idxs;
    if (exist("breakpoint_mask") == 1) % i.e. is an existing variable
      model_idxs &= breakpoint_mask;
    end % if
    fit           = polyfit(t(model_idxs), data.timeseries(model_idxs), 1);
    ms(i)         = fit(1);
    cs(i)         = fit(2);
    y_orig        = data.timeseries(flatten_idxs);
    x             = t(flatten_idxs);
    analysed.flattened(flatten_idxs) ...
                  = y_orig - (ms(i)*x + cs(i));
  end % for

  analysed.ms = ms;
  analysed.cs = cs;

end % function
