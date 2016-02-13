function plot_profile()

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  plot_name = "profile";

  a = create_axes(plot_name);
  if (isempty(a))
    return
  end % if
  ax = plots.(plot_name).axis;

  % Check that period has been set
  if (isempty(analysis.period))
    errordlg("The period has not been set");
    delete(figures.(plot_name).fig_handle);
    figures.(plot_name).fig_handle = [];
    figures.(plot_name).ax_handle = [];
    return
  end % if

  % Check that number of profile bins has been set
  if (isempty(analysis.nprofile_bins))
    errordlg("The number of profile bins has not been set");
    delete(figures.(plot_name).fig_handle);
    figures.(plot_name).fig_handle = [];
    figures.(plot_name).ax_handle = [];
    return
  end % if

  % Calculate profile
  calc_profile();

  % Plot it up!
  phase = [0:(analysis.nprofile_bins - 1)] / analysis.nprofile_bins;
  plot(a, phase, analysed.profile);
  xlabel(a, "Phase");
  ylabel(a, "Flux (UNITS)");
  set_title(plot_name);

  % Get/Set axis limits
  if (~isempty(ax))
    axis(a, ax);
  else
    plots.(plot_name).axis = axis(a);
  end % if

  % Shade the masked area
  if (~isempty(analysis.profile_mask))
    hold on;
    x1 = analysis.profile_mask(1);
    x2 = analysis.profile_mask(2);
    ymin = plots.(plot_name).axis(3);
    ymax = plots.(plot_name).axis(4);
    mygreen = [0.5, 1.0, 0.5];
    if (x1 <= x2)
      area(a, analysis.profile_mask, [ymax,ymax], ymin, "facecolor", mygreen, "linewidth", 0);
    else
      area(a, [0,x2], [ymax,ymax], ymin, "facecolor", mygreen, "linewidth", 0);
      area(a, [x1,1], [ymax,ymax], ymin, "facecolor", mygreen, "linewidth", 0);
    end % if
    hold off;
  end % if

end % function
