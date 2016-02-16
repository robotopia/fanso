function plot_tdfs()

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  plot_name = "tdfs";

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

  % Stack the pulses and perform 2D FFT
  stackpulses(true); % <-- "true" = do only the visible timeseries
  analysed.tdfs.original = fft2(analysed.stacked.transformed);

  % Shift in x and y so that (shifted) DC is in centre of grid (as opposed to position (1,1))
  s = size(analysed.tdfs.original);
  nxs    = s(2);                 nys    = s(1);
  xmin   = -floor((nxs-1)/2);    xmax   = ceil((nxs-1)/2);
  ymin   = -floor((nys-1)/2);    ymax   = ceil((nys-1)/2);
  xshift = -xmin;                yshift = -ymin;

  analysed.tdfs.xs      = [xmin:xmax];          % Units of   v_l*P1/(2*pi)   or   v_l*P1?
  analysed.tdfs.ys      = [ymin:ymax] / nys;    % Units v_t * P1
  ymin ./= nys;
  ymax ./= nys;

  analysed.tdfs.centred = shift(analysed.tdfs.original, xshift, 2);
  analysed.tdfs.centred = shift(analysed.tdfs.centred,  yshift, 1);

  % Apply filters
  apply_tdfs_filters();

  to_be_plotted = abs(analysed.tdfs.filtered);
  clabel_text   = 'Amplitude';

  % Are we plotting power instead?
  if (plots.(plot_name).ispower)
    to_be_plotted .^= 2;
    clabel_text = 'Power';
  end % if

  % Is it log?
  if (plots.(plot_name).islog)
    to_be_plotted = log10(to_be_plotted);
    clabel_text = ["log_{10} (", clabel_text, ")"];
  end % if

  % Plot it up!
  plots.(plot_name).autoscale = [xmin, xmax, ymin, ymax];
  imagesc(a, analysed.tdfs.xs, analysed.tdfs.ys, to_be_plotted);
  axis("xy");

  apply_colormap(plot_name);

  xlabel(a, "2*pi*v_1*P_1");
  ylabel(a, "v_1P_1");
  set_title(plot_name);

  % If in filter edit mode, plot filter boundaries,
  % otherwise plot shifted DC axes
  if (~isfield(analysed, "filter")) % = not in filter edit mode
    % Plot the point of the new centre
    if (~isempty(analysis.shift_DC))
      hold(a, "on");
      X = [analysis.shift_DC(1), xmin;
           analysis.shift_DC(1), xmax];
      Y = [ymin, analysis.shift_DC(2);
           ymax, analysis.shift_DC(2)];
      plot(a, X, Y, "g");
      hold(a, "off");
    end % if

  else % = in filter edit mode

    % Plot the filters if in filter mode
    hold(a, "on");
    for n = 1:rows(analysis.filters) % For each filter in the list...

      centre    = analysis.filters(n,1);
      width     = analysis.filters(n,2);
      direction = analysis.filters(n,3);

      switch direction

        case 0 % Horizontal filter

          Xouter = [xmin, xmin;
                    xmax, xmax];
          Xinner = Xouter;
          Youter = [centre+width,   centre-width;
                    centre+width,   centre-width];
          Yinner = [centre+width/2, centre-width/2;
                    centre+width/2, centre-width/2];

        case 1 % Vertical filter

          Xouter = [centre+width,   centre-width;
                    centre+width,   centre-width];
          Xinner = [centre+width/2, centre-width/2;
                    centre+width/2, centre-width/2];
          Youter = [ymin, ymin;
                    ymax, ymax];
          Yinner = Youter;

      end % switch

      plot(a, Xouter, Youter, "g--", Xinner, Yinner, "g");

    end % for
    hold(a, "off");
  end % if

  % Add a colorbar
  c = colorbar(a, "ylabel", clabel_text);
  %colorbar(c, "ylabel", clabel_text);

  % (For some reason, imagesc() wants to change the axes position, so we don't let it)
  set(a, "position", [0.13, 0.11, 0.62, 0.815]);
  set(c, "position", [0.8275, 0.11, 0.0465, 0.815]);

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

function apply_tdfs_filters()

  global analysis;
  global analysed;

  [X, Y] = meshgrid(analysed.tdfs.xs, analysed.tdfs.ys);

  % --note-- I'm not sure this algorithm covers edge cases appropriately
  nfilters = rows(analysis.filters);
  transmission_total = ones(size(analysed.tdfs.centred));
  for n = 1:nfilters
    switch analysis.filters(n,3)
      case 0 % horizontal filter
        transmission = 2*abs(Y - analysis.filters(n,1))/analysis.filters(n,2) - 1;
      case 1 % vertical filter
        transmission = 2*abs(X - analysis.filters(n,1))/analysis.filters(n,2) - 1;
    end % switch
    transmission(transmission > 1) = 1;
    transmission(transmission < 0) = 0;
    transmission_total = min(transmission_total, transmission);
  end % for

  % Multiply the filter by the 2DFS
  analysed.tdfs.filtered = analysed.tdfs.centred .*= transmission_total;

end % function


