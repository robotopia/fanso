function plot_fft()

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  plot_name = "fft";

  a = create_axes(plot_name);
  if (isempty(a))
    return
  end % if
  ax = plots.(plot_name).axis;

  % Apply windowing, etc.
  apply_transforms();

  % Calculate values for the FFT abcissa...
  n                        = length(analysed.transformed);
  df                       = 1/(analysed.dt*n);
  analysed.spectrum_freqs  = [0:(n-1)] * df;
  % ...and the FFT ordinate
  analysed.spectrum_vals   = fft(analysed.transformed);

  to_be_plotted = abs(analysed.spectrum_vals);
  ylabel_text   = "Amplitude";

  % Are we plotting power instead?
  if (plots.(plot_name).ispower)
    to_be_plotted .^= 2;
    ylabel_text = "Power";
  end % if

  % Set plot function according to whether log plots are requested
  if (plots.(plot_name).islog)
    plot_fcn = @semilogy;
  else
    plot_fcn = @plot;
  end % if

  % Plot up the FFT
  plot_fcn(a, analysed.spectrum_freqs, to_be_plotted, 'b');
  xlabel(a, ["Frequency (", data.frequnits, ")"]);
  ylabel(a, ylabel_text);
  set_title(plot_name);

  % Get/Set axis limits
  if (~isempty(ax))
    axis(a, ax);
  else
    plots.(plot_name).axis = axis(a);
  end % if

end % function
