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

  % Calculate the FFT
  calc_fft();
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

  % Create an informative title
  if (isfield(analysed, "fft_title"))
    fft_title = analysed.fft_title;
  else
    if (analysis.only_visible)
      tax = plots.timeseries.axis;
      xlo = num2str(max([tax(1), 0]));
      xhi = num2str(min([tax(2), analysed.t(end)]));
    else
      xlo = "0";
      xhi = num2str(analysed.N * analysed.dt);
    end % if
    fft_title = ["FFT of timeseries between ", xlo, " and ", xhi, " ", data.timeunits];
    if (analysis.apply_hamming)
      fft_title = [fft_title, "\nHamming window applied"];
    end % if
    if (analysis.apply_hanning)
      fft_title = [fft_title, "\nHanning window applied"];
    end % if
    if (analysis.zeromean)
      fft_title = [fft_title, "\nZero-meaning timeseries"];
    end % if
    if (analysis.zeropad)
      fft_title = [fft_title, "\nZero-padding timeseries"];
    end % if
  end % if
  title(a, fft_title);

  % Get/Set axis limits
  if (~isempty(ax))
    axis(a, ax);
  else
    plots.(plot_name).axis = axis(a);
  end % if

end % function
