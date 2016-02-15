function plot_pulsestack()

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  plot_name = "pulsestack";

  a = create_axes(plot_name);
  if (isempty(a))
    return
  end % if
  ax = plots.(plot_name).axis;
  cax = plots.(plot_name).cax;

  % Make sure the period has been set
  if (isempty(analysis.period))
    errordlg("The period has not been set");
    delete(figures.(plot_name).fig_handle);
    figures.(plot_name).fig_handle = [];
    figures.(plot_name).ax_handle  = [];
    return
  end % if

  % Stack the pulses!
  stackpulses(false); % <-- "false" = do the whole timeseries

  % Get appropriate values for x and y axes
  nxs = columns(analysed.stacked.all);
  nys =    rows(analysed.stacked.all);

  xs = [0:(nxs-1)] / nxs;  % Phase
  ys = [1:nys];            % Pulse number

  % Plot it up!
  plots.(plot_name).autoscale = [0, 1, 0.5, nys+0.5];
  imagesc(a, xs, ys, analysed.stacked.all);
  apply_colormap(plot_name);
  axis("xy");

  % Draw horizontal lines if only_visible is on
  if (analysis.only_visible)
    ax = plots.timeseries.axis;
    ymin = max([ax(1)/analysis.period,0.5]);
    ymax = min([ax(2)/analysis.period,nys+0.5]);
    X = [-0.5,nxs-0.5] / nxs;
    Y1 = [ymin, ymin];
    Y2 = [ymax, ymax];
    hold on;
    plot(a, X,Y1,'g',X,Y2,'g');
    hold off;
  end % if

  % Add a colorbar
  colorbar(a);

  xlabel(a, "Phase");
  ylabel(a, "Pulse number");
  set_title(plot_name);

  % Get/Set axis limits
  if (~isempty(ax))
    axis(a, ax);
  else
    plots.(plot_name).axis = axis(a);
  end % if

  % Get/Set caxis limits
  if (~isempty(cax))
    caxis(a, cax);
  else
    plots.(plot_name).cax  = caxis(a);
  end % if

end % function
