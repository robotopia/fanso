function flatten()

  global data;
  global analysis;
  global analysed;

  if (isempty(analysis.breakpoints))
    return
  end % if

  % Set up flattened vector
  analysed.flattened = zeros(size(data.timeseries));

  % Book-end the breakpoints with initial and final values
  N = length(data.timeseries);
  dt = 1/data.samplingrate;
  t = [0:(N-1)]' * dt;

  bps = [t(1)-1, analysis.breakpoints, t(end)+1];

  % Loop through the breakpoints and detrend each segment
  ms = zeros(length(bps)-1, 1);
  cs = zeros(length(bps)-1, 1);
  for i = 1:(length(bps)-1)
    flatten_idxs  = (t >= bps(i)) & (t < bps(i+1));
    model_idxs    = flatten_idxs;
    if (isfield(analysed, "breakpoint_mask"))
      model_idxs &= analysed.breakpoint_mask;
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

