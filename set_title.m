function set_title(plot_name)

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  a = figures.(plot_name).ax_handle;

  if (~isempty(a))
    analysed_title = sprintf("%s_title", plot_name);
    if (isfield(analysed, analysed_title))
      title(a, analysed.(analysed_title));
    else
      switch plot_name
        case "timeseries"
          title(a, sprintf("Sampling rate = %f", data.samplingrate));

        case "fft"
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
          title(a, fft_title);

        case "profile"
          title(a, sprintf("Period = %.12f %s;    No. of bins = %d", analysis.period, data.timeunits, analysis.nprofile_bins));

        case "pulsestack"
          title(a, "");

        case "tdfs"
          title(a, "");

        case "modenv"
          title(a, "");

      end % switch
    end % if
  end % if

end % function
