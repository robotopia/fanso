function plot_timeseries()

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  plot_name = "timeseries";

  a = create_axes(plot_name);
  if (isempty(a))
    return
  end % if
  ax = plots.(plot_name).axis;

  % "analysed" structure (to avoid recalculating the same things more than once)
  if (isempty(analysed))
    analysed.N  = length(data.timeseries);             % The length of the timeseries
    analysed.dt = 1/data.samplingrate;                 % The time between adjacent samples
    analysed.t  = [0:(analysed.N-1)]' * analysed.dt;   % The time axis
  end % if

  % Plot different things depending on whether the breakpoints are
  % to be "applied" or not (i.e. whether the timeseries has been flattened).
  if (~isempty(analysis.breakpoints))
    flatten();
  end % if

  if (analysis.apply_bps)
    plot(a, analysed.t, analysed.flattened, 'b');
  else
    if(~isempty(analysis.breakpoints))

      % Prepare vertical lines for breakpoint plotting
      bp_xs = [1;1] * analysis.breakpoints;
      ymin = min(data.timeseries);
      ymax = max(data.timeseries);
      bp_ys = repmat([ymin;ymax],size(analysis.breakpoints));

      % Prepare the detrend lines for plotting
      model_xs = [0, analysis.breakpoints;
                  analysis.breakpoints, analysed.t(end)];
      m_mat = repmat(analysed.ms',2,1);
      c_mat = repmat(analysed.cs',2,1);
      model_ys = model_xs .* m_mat + c_mat;

      % Plot the original timeseries, breakpoints, and detrend lines
      plot(a, analysed.t, data.timeseries, 'b', ...
              bp_xs,      bp_ys,           'r', "linewidth", 2.0, ...
              model_xs,   model_ys,        'g', "linewidth", 2.0);
    else
      % Plot the timeseries!
      plot(a, analysed.t, data.timeseries, 'b');
    end % if
  end % if

  xlabel(a, ["Time (", data.timeunits, ")"]);
  ylabel(a, "Flux (UNITS)"); % <-- add a field for this too

  % Get/Set axis limits
  if (~isempty(ax))
    axis(a, ax);
  else
    plots.(plot_name).axis = axis(a);
  end % if

end % function
