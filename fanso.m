function fanso(init_filename)
% -- Function: fanso ()
%     Analyse a timeseries using the Fourier ANalysis Suite for Octave.
%     Requires FLTK.
%
%     Written by Sam McSweeney, 2016, Creative Commons Licence
%     sammy.mcsweeney@gmail.com

  % First, check if there is already an instance of FANSO open
  global __FANSO_INSTANCE__
  if (isempty(__FANSO_INSTANCE__)) % i.e. there is NO instance already open
    __FANSO_INSTANCE__ = 1;
  else
    % Exit gracefully
    disp('FANSO already running');
    return
  end

  % Require that we are using the correct graphics toolkit (FLTK)
  gtk = 'fltk';
  if (~any(strcmp(available_graphics_toolkits(), gtk)))
    disp(['This function requires ',gtk,' to be installed.'])
    return;
  else
    graphics_toolkit(gtk);
  end % if

  % No file is loaded yet...
  global filename
  filename = [];

  % Setup figure and plot handling system
  global figures;
  figures = struct();

  plot_names = {"timeseries", "fft", "profile", "hrfs", "waterfall", "tdfs", "modenv"};
  plot_dims  = [1, 1, 1, 2, 2, 2, 1];
  plot_isfft = [0, 1, 0, 1, 0, 1, 0];

  for n = 1:length(plot_names)
    figures.(plot_names{n}).drawfcn     = str2func(["plot_", plot_names{n}]); % Function that does the plotting
    figures.(plot_names{n}).fig_handle  = [];                                 % Handle to the figure
    figures.(plot_names{n}).ax_handle   = [];                                 % Handle(s) to the plot axes
    figures.(plot_names{n}).quantise    = 0;                                  % Whether a click is rounded to the nearest "image pixel"
    figures.(plot_names{n}).dims        = plot_dims(n);                       % Number of dimensions in plot
    figures.(plot_names{n}).isfft       = plot_isfft(n);                      % Whether the plot is of type "FFT"
  end % for

  % Set the default positions for the timeseries and fft figures
  screensize = get(0, 'screensize');
  gapsize    = 100;

  winsize_x  = screensize(3) - 2*gapsize;
  winsize_y  = floor((screensize(4) - 3*gapsize)/2);
  winpos_x   = gapsize;
  winpos_y   = 2*gapsize + winsize_y;
  figures.timeseries.defaultpos = [winpos_x, winpos_y, winsize_x, winsize_y];

  winpos_y   = gapsize;
  figures.fft.defaultpos = [winpos_x, winpos_y, winsize_x, winsize_y];

  % Load known pulsar periods
  load_periods();

  % Draw the main (=timeseries) window
  create_figure("timeseries");

  % Initially, ensure that there are no changes to be saved
  set_unsaved_changes(false);

  % Load initial file, if given
  if (nargin == 1)
    load_fan(init_filename);
  end % if

end % function


function load_periods()

  global pulsars

  f = fopen('periods.dat');

  % First, count the lines
  nlines = fskipl(f, Inf);
  frewind(f);

  % Now read them in
  pulsars = struct();
  for n = 1:nlines
    [pulsarname, pulsarperiod] = fscanf(f, '%s %f', "C");
    pulsars.(pulsarname).period = pulsarperiod;
  end % for

end % function

%%%%%%%%%%%%%%%%%%%%%
% Everything after this point will eventually be relocated
%

function Xtoggle_flatten(src, data)

  global fig_handles; % <-- is this (and other similar lines) needed?

  global apply_bps
  apply_bps = ~apply_bps;

  set_unsaved_changes(true);

  if (apply_bps)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  get_axes();
  replot();

end % function

function Xcreate_figureX(src, data, fig_no)

  global fig_handles

  switch fig_no

    case 2 % The FFT figure
      global screensize
      global gapsize

      winsize_x  = screensize(3) - 2*gapsize;
      winsize_y  = floor((screensize(4) - 3*gapsize)/2);
      winpos_x   = gapsize;
      winpos_y   = gapsize;

      fig_handles(fig_no) = figure("Name", "Fourier Transform", ...
                                   "Position", [winpos_x, winpos_y, winsize_x, winsize_y], ...
                                   "CloseRequestFcn", {@close_figure, fig_no});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up menu for  FFT figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      m_plot           = create_plot_menu(fig_no);

      m_analyse        = uimenu('label', '&Analyse');
      m_analyse_period = uimenu(m_analyse, 'label', '&Select period (P1)', 'callback', @select_period);

    case 3 % The profile figure
      fig_handles(fig_no) = figure("Name", "Profile", ...
                                   "CloseRequestFcn", {@close_figure, fig_no});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up menu for profile figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      create_plot_menu(fig_no);

      m_mask        = uimenu('label', '&Mask');
      m_mask_clear  = uimenu(m_mask, 'label', '&Clear', 'callback', @clear_mask);
      m_mask_select = uimenu(m_mask, 'label', '&Select', 'callback', @select_mask);

    case 4 % The harmonic resolved fluctuation spectrum figure
      fig_handles(fig_no) = figure("Name", "Harmonic Resolved Fluctuation Spectrum", ...
                                   "CloseRequestFcn", {@close_figure, fig_no});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up menu for HRFS figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      create_plot_menu(fig_no);
      create_colormap_menu(fig_no);

    case 5 % The waterfall plot of the timeseries
      fig_handles(fig_no) = figure("Name", "Waterfall plot", ...
                                   "CloseRequestFcn", {@close_figure, fig_no});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up menu for waterfall figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      global m_waterfall_colormap
      m_plot     = create_plot_menu(fig_no);
      m_plot_2D  = uimenu(m_plot, 'label', '2D', 'accelerator', '2', 'separator', 'on', ...
                                  'callback', {@set_waterfall_plot_type, 0});
      m_plot_3D  = uimenu(m_plot, 'label', '3D', 'accelerator', '3', ...
                                  'callback', {@set_waterfall_plot_type, 1});
      m_waterfall_colormap = create_colormap_menu(fig_no);

    case 6 % The 2DFS plot
      fig_handles(fig_no) = figure("Name", "2D Fluctuation Spectrum", ...
                                   "CloseRequestFcn", {@close_figure, fig_no});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up menu for 2DFS figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      create_plot_menu(fig_no);
      create_colormap_menu(fig_no);
      m_tdfs_analyse         = uimenu('label', 'Analyse');
      m_tdfs_analyse_peaks   = uimenu(m_tdfs_analyse, 'label', 'Show peaks', 'callback', @toggle_peaks);
      m_tdfs_analyse_p2p3    = uimenu(m_tdfs_analyse, 'label', 'Choose P2, P3', 'callback', {@click_p2p3, fig_no});
      m_tdfs_analyse_hfilt   = uimenu(m_tdfs_analyse, 'label', 'Create horizontal filter', 'separator', 'on', ...
                                      'callback', {@create_filter, fig_no, "horizontal"});
      m_tdfs_analyse_vfilt   = uimenu(m_tdfs_analyse, 'label', 'Create vertical filter', ...
                                      'callback', {@create_filter, fig_no, "vertical"});
      m_tdfs_analyse_delfilt = uimenu(m_tdfs_analyse, 'label', 'Delete filter', 'separator', 'on', ...
                                      'callback', {@delete_filter, fig_no});
      m_tdfs_analyse_delall  = uimenu(m_tdfs_analyse, 'label', 'Delete all filters', 'enable', 'off');
      m_tdfs_analyse_shiftDC = uimenu(m_tdfs_analyse, 'label', 'Shift DC', 'separator', 'on', ...
                                      'callback', {@click_shiftDC, fig_no});
      m_tdfs_analyse_modenv = uimenu(m_tdfs_analyse, 'label', 'Plot complex envelopes', 'separator', 'on', 'callback', @plot_modenv);

    case 7 % The complex envelope = the result of E&S (2002) analysis
      fig_handles(fig_no) = figure("Name", "Complex Envelope", ...
                                   "CloseRequestFcn", {@close_figure, fig_no});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up menu for Complex Envelope figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Empty menu %

  end % switch

end % function

function Xm_plot = create_plot_menu(fig_no, parent = 0)

  global fig_handles
  figure(fig_handles(fig_no));

  % Plot menu
  if (parent)
    m_plot = uimenu(parent, 'label', 'Plot');
  else
    m_plot = uimenu('label', 'Plot');
  end % if

  % Children menu items
  m_plot_save = uimenu(m_plot, 'label', 'Redraw all plots', 'accelerator', '1', 'callback', {@update_plots, fig_no});

  % For FFT plots, allow option of plotting absolute values or power
  if (any(fig_no == [2,4,6])) % 2, 4, 6 are the FFT plots
    m_plot_abs   = uimenu(m_plot, 'label', 'Amplitudes', 'separator', 'on', 'callback', {@set_fft_plot_type, fig_no, nan});
    m_plot_power = uimenu(m_plot, 'label', 'Power', 'callback', {@set_fft_plot_type, fig_no, 2});
  end

  % Allow the option to turn on/off log-plotting for select axes depending on the plot type
  switch fig_no
    case {2} % FFT plot
      m_plot_xlog = uimenu(m_plot, 'label', 'Log x-axis', 'separator', 'on', 'callback', {@toggle_log, fig_no, 'x'});
      m_plot_ylog = uimenu(m_plot, 'label', 'Log y-axis', 'callback', {@toggle_log, fig_no, 'y'});
    case {4,6} % HRFS and 2DFS plots
      m_plot_clog = uimenu(m_plot, 'label', 'Log pixel values', 'separator', 'on', 'callback', {@toggle_log, fig_no, 'c'});
  end % switch

end % function

function Xm_colormap = create_colormap_menu(fig_no, plot_no)

  global plot_params
  global fig_handles
  figure(fig_handles(fig_no));

  % Parent menu
  m_colormap = uimenu('label', '&Colormap');

  % Children menu items
  m_type = uimenu(m_colormap, 'label', 'Type');

  type_list = colormap("list");
  for n = 1:length(type_list)
    uimenu(m_type, 'label', type_list{n}, 'callback', {@set_colormap_type, fig_no, n});
  end % for

  uimenu(m_colormap, 'label', 'Invert', 'checked', on_off(plot_params(fig_no, 18)), ...
                     'callback', {@toggle_plot_param, plot_no, "cinv"});
  uimenu(m_colormap, 'label', 'Set dynamic range limits', 'separator', 'on', ...
                     'callback', {@scale_caxis, plot_no, 0, 0});
  uimenu(m_colormap, 'label', 'Increase Max by 10% of range', 'accelerator', '+', ...
                     'callback', {@scale_caxis, plot_no, 0,  0.1});
  uimenu(m_colormap, 'label', 'Decrease Max by 10% of range', 'accelerator', '_', ...
                     'callback', {@scale_caxis, plot_no, 0, -0.1});
  uimenu(m_colormap, 'label', 'Increase Min by 10% of range', 'accelerator', '=', ...
                     'callback', {@scale_caxis, plot_no,  0.1, 0});
  uimenu(m_colormap, 'label', 'Decrease Min by 10% of range', 'accelerator', '-', ...
                     'callback', {@scale_caxis, plot_no, -0.1, 0});

end % function

function Xupdate_plots(src, data)

  get_axes();
  replot();

end % function

function Xset_waterfall_plot_type(src, data, newvalue)

  global waterfall_plot_type
  global m_waterfall_colormap

  if (waterfall_plot_type ~= newvalue)

    % Change the value!
    waterfall_plot_type = newvalue;
    set_unsaved_changes(true);

    % Enable/Disable the colormap menu for 2D/3D waterfall plots
    set(m_waterfall_colormap, "enable", on_off(~newvalue)); % 0 = 2D = menu ON;  1 = 3D = menu OFF

    % Update the figure
    get_axes(0,0,5);
    replot(5);

  end % if

end % function

function Xset_axis(fig_no, axis_char, newmin = NaN, newmax = NaN)

  global plot_params

  try
    colmax = 2*find(axis_char == "xyzc");
    colmin = colmax - 1;
  catch
    errordlg({["Unfamiliar axis code (\"", axis_char, "\":"], lasterr()});
    return
  end % try_catch

  % Get current values
  curr_min = plot_params(fig_no, colmin);
  curr_max = plot_params(fig_no, colmax);

  % If new values are NaNs, then ask for values via a dialog box
  switch isnan([newmin, newmax])

    case [true, true]

      % Ask for new values via a dialog box
      cstrs = inputdlg({"Min:", "Max:"}, ["Enter new ", axis_char, "-axis limits"], ...
                       1, {num2str(curr_min), num2str(curr_max)});

      if (isempty(cstrs)) % They pushed cancel
        return
      end % if

      newmin = str2num(cstrs{1});
      newmax = str2num(cstrs{2});

    case [true, false]

      % Ask for one new value
      cstr = inputdlg({"Min:"}, ["Enter new ", axis_char, "-axis minimum"], ...
                      1, {num2str(curr_min)});

      if (isempty(cstr)) % They pushed cancel
        return
      end % if

      newmin = str2num(cstr{1});

    case [false, true]

      % Ask for one new value
      cstr = inputdlg({"Max:"}, ["Enter new ", axis_char, "-axis maximum"], ...
                      1, {num2str(curr_max)});

      if (isempty(cstr)) % They pushed cancel
        return
      end % if

      newmax = str2num(cstr{1});

  end % switch

  % Check if new values are different from old
  if ((newmin ~= curr_min) || (newmax ~= curr_max))
    set_unsaved_changes(true);
  end % if

  % Actually make the changes
  try
    plot_params(fig_no, colmin) = newmin;
    plot_params(fig_no, colmax) = newmax;
  catch
    errordlg({"Only numeric values allowed:", lasterr()});
    return
  end % try_catch

  replot(fig_no);

end % function

function Xscale_caxis(src, data, fig_no, minfactor, maxfactor)

  if (all(~[minfactor, maxfactor])) % Putting in zeros for these parameters brings up a dialog box

    set_axis(fig_no, 'c');

  else

    % Set new values as scaled version of the old ones
    global plot_params

    curr_cmin = plot_params(fig_no, 7);
    curr_cmax = plot_params(fig_no, 8);

    crange = curr_cmax - curr_cmin;

    newmin = curr_cmin + crange*minfactor;
    newmax = curr_cmax + crange*maxfactor;

    set_axis(fig_no, 'c', newmin, newmax);

  end % if

end % function

function Xstr = on_off(bool)
% Converts 1,0 to "on","off" respectively
  if (bool)
    str = "on";
  else
    str = "off";
  end % if
end % function

function Xplot_fft(src, data)

  global fig_handles

  global spectrum_freqs
  global spectrum_vals

  global plot_params

  % Switch to/Create FFT figure and keep track of the view window
  fig_no = 2;
  first_time = (fig_handles(fig_no) == 0);

  if (first_time)
    create_figure(0,0,fig_no);
  end % if

  figure(fig_handles(fig_no));

  % Calculate the FFT
  calc_fft();
  to_be_plotted = abs(spectrum_vals);
  ylabel_text   = 'Amplitude';

  % Are we plotting power instead?
  if (plot_params(fig_no, 14) == 2)
    to_be_plotted .^= 2;
    ylabel_text = 'Power';
  end % if

  % Set plot function according to whether log plots are requested
  switch [plot_params(fig_no, 9:10)]
    case [0,0]
      plot_fcn = @plot;
    case [1,0]
      plot_fcn = @semilogx;
    case [0,1]
      plot_fcn = @semilogy;
    case [1,1]
      plot_fcn = @loglog;
    otherwise
      errordlg("Non-logical values for plot log type");
      return
  end % switch

  % Plot up the FFT
  figure(fig_handles(fig_no));
  plot_fcn(spectrum_freqs, to_be_plotted, 'b');
  xlabel('Frequency (Hz)');
  ylabel(ylabel_text);

  % Get/Set axis limits
  if (~set_axes(fig_no))
    get_axes(0,0,fig_no);
  end % if

end % function

function Xselect_period(src, data)

  global fig_handles
  figure(fig_handles(2));

  global period;

  title("Click on a harmonic of the desired frequency");

  [x, y, button] = ginput(1);

  % Get the harmonic number via a dialog box
  cstr = inputdlg({"Harmonic number of selected point:"}, "Harmonic", 1, {"1"});
  if (~isempty(cstr))
    nharm = str2num(cstr{1});
    newperiod = nharm / x;
    set_period(newperiod);
  end % if

  % Reset title
  figure(fig_handles(2));
  title("");

end % function

function Xplot_profile(src, data)

  global fig_handles

  global profile
  global nprofile_bins
  global period
  global profile_mask

  if (isempty(period))
    errordlg('The period has not been set');
    return
  end % if

  if (isempty(nprofile_bins))
    errordlg('The number of profile bins has not been set');
    return
  end % if

  % Switch to/Create FFT figure and keep track of the view window
  fig_no = 3;
  first_time = (fig_handles(fig_no) == 0);

  if (first_time)
    create_figure(0,0,fig_no);
  end % if

  h = fig_handles(fig_no);

  % Calculate profile
  calc_profile();

  % Plot it up!
  figure(h);
  hold off;
  dphase = 1/nprofile_bins;
  phase = [0:dphase:(1-dphase/2)]; % Ensure that there are the correct number of bins
  plot(phase, profile);
  xlabel('Phase');
  ylabel('Flux (arbitrary units)');
  profile_title = sprintf('Period = %.12f s;    No. of bins = %d', period, nprofile_bins);
  title(profile_title);

  % Get/Set axis limits
  if (~set_axes(fig_no))
    get_axes(0,0,fig_no);
  end % if

  % Shade the masked area
  figure(h);
  ax = axis();
  ymin = ax(3);
  ymax = ax(4);
  hold on;
  mygreen = [0.5,1.0,0.5];
  if (profile_mask(1) <= profile_mask(2))
    area(profile_mask, [ymax, ymax], ymin, "FaceColor", mygreen, "linewidth", 0);
  else
    area([0,profile_mask(2)], [ymax,ymax], ymin, "FaceColor", mygreen, "linewidth", 0);
    area([profile_mask(1),1], [ymax,ymax], ymin, "FaceColor", mygreen, "linewidth", 0);
  end % if

end % function

function Xplot_hrfs(src, data)

  global fig_handles

  global dt
  global period

  global zeropad

  global spectrum_freqs
  global spectrum_vals

  global plot_params

  % Switch to/Create HRFS figure and keep track of the view window
  fig_no = 4;
  first_time = (fig_handles(fig_no) == 0);

  if (first_time)
    create_figure(0,0,fig_no);
  end % if

  figure(fig_handles(fig_no));

  % Turn on zero-padding
  if (~zeropad)
    toggle_zeropad();
  end % if

  % Calculate the FFT
  calc_fft();
  to_be_plotted = abs(spectrum_vals);
  clabel_text   = 'Amplitude';

  % Are we plotting power instead?
  if (plot_params(fig_no, 16) == 2)
    to_be_plotted .^= 2;
    clabel_text = 'Power';
  end % if

  % Prepare the matrix of values
  N             = length(to_be_plotted);
  nx            = round(N*dt/period);
  ny            = floor(N/(2*nx));
  n             = nx*ny;
  to_be_plotted = to_be_plotted(1:n);
  to_be_plotted = reshape(to_be_plotted,nx,ny)';

  xs = [0:(nx-1)]*spectrum_freqs(2)*period;
  ys = [0:(ny-1)];

  % Plot log values, if requested
  if (plot_params(fig_no, 12))
    to_be_plotted = log10(to_be_plotted);
    clabel_text = ['log_{10} (', clabel_text, ')'];
  end % if

  % Plot up the FFT
  figure(fig_handles(fig_no));
  imagesc(xs, ys, to_be_plotted);
  xlabel('Frequency*P');
  ylabel('Harmonic Number');
  axis("xy");
  colorbar('ylabel', clabel_text);
  apply_colormap(fig_no);

  % Get/Set axis limits
  if (~set_axes(fig_no))
    get_axes(0,0,fig_no);
  end % if

end % function

function Xplot_waterfall(src, data)

  global fig_handles

  global period
  global timeseries_grid
  global plot_params

  global waterfall_plot_type
  global only_visible

  % Make sure the period has been set
  if (isempty(period))
    errordlg('The period has not been set');
    return
  end % if

  % Switch to/Create HRFS figure and keep track of the view window
  fig_no = 5;
  first_time = (fig_handles(fig_no) == 0);

  if (first_time)
    create_figure(0,0,fig_no);
  end % if

  figure(fig_handles(fig_no));

  reshape_timeseries_into_grid(false); % <-- false = do whole timeseries, not just visible part

  % Get appropriate values for x and y axes
  nxs = columns(timeseries_grid);
  nys =    rows(timeseries_grid);

  xs = [0:(nxs-1)] / nxs;  % Phase
  ys = [1:nys];            % Pulse number

  switch waterfall_plot_type
    case 0 % 2D
      imagesc(xs, ys, timeseries_grid);
      apply_colormap(fig_no);
      axis("xy");

      if (only_visible)
        ymin = max([plot_params(1,1)/period,0.5]);
        ymax = min([plot_params(1,2)/period,nys+0.5]);
        hold on;
        X = [-0.5,nxs-0.5];
        Y1 = [ymin, ymin];
        Y2 = [ymax, ymax];
        plot(X,Y1,'g',X,Y2,'g');
        hold off;
      end % if

      colorbar();

    case 1 % 3D
      [Xs, Ys] = meshgrid(xs, ys);
      waterfall(Xs, Ys, timeseries_grid);
      colormap([0,0,0]); % i.e. all black lines
    otherwise
      error("Unknown waterfall plot type requested");
  end % switch

  xlabel('Phase');
  ylabel('Pulse number');

  % Get/Set axis limits
  if (~set_axes(fig_no))
    get_axes(0,0,fig_no);
  end % if

end % function

function Xplot_tdfs(src, data)

  global fig_handles
  global plot_params

  global timeseries_grid
  global tdfs

  global period
  global P2hat
  global P3hat
  global show_peaks
  global filters
  global new_filter
  global shift_DC

  % Switch to/Create 2DFS figure and keep track of the view window
  fig_no = 6;
  first_time = (fig_handles(fig_no) == 0);

  if (first_time)
    create_figure(0,0,fig_no);
  end % if

  figure(fig_handles(fig_no));

  % Calculate the 2DFS
  reshape_timeseries_into_grid();
  tdfs = fft2(timeseries_grid);

  apply_tdfs_filters();
  apply_shift_DC();

  to_be_plotted = abs(tdfs);
  clabel_text   = 'Amplitude';

  % Are we plotting power instead?
  if (plot_params(fig_no, 16) == 2)
    to_be_plotted .^= 2;
    clabel_text = 'Power';
  end % if

  % Shift in x and y so that (shifted) DC is in centre of grid (as opposed to position (1,1)
  nxs = columns(to_be_plotted);  nys = rows(to_be_plotted);
  xmin = -floor((nxs-1)/2);      xmax = ceil((nxs-1)/2);
  ymin = -floor((nys-1)/2);      ymax = ceil((nys-1)/2);
  xshift = -xmin;                yshift = -ymin;

  to_be_plotted = shift(to_be_plotted, xshift, 2);
  to_be_plotted = shift(to_be_plotted, yshift);

  xs = [xmin:xmax];          % Units of   v_l*P1/(2*pi)   or   v_l*P1?
  ys = [ymin:ymax] / nys;    % Units v_t * P1

  % Plot log values, if requested
  if (plot_params(fig_no, 12))
    to_be_plotted = log10(to_be_plotted);
    clabel_text = ['log_{10} (', clabel_text, ')'];
  end % if

  % Draw plot
  imagesc(xs, ys, to_be_plotted);
  axis("xy");

  % Show local maxima
  if (show_peaks)
    % Calculate local maxima positions
    m1 = to_be_plotted >= shift(to_be_plotted,  1, 1); % Pixel is >= the one below it
    m2 = to_be_plotted >= shift(to_be_plotted, -1, 1); % Pixel is >= the one above it
    m3 = to_be_plotted >= shift(to_be_plotted,  1, 2); % Pixel is >= the one to the right of it
    m4 = to_be_plotted >= shift(to_be_plotted, -1, 2); % Pixel is >= the one to the left of it

    % Apply threshold so that we don't get every single peak in the noise
    threshold = 0.10; % i.e. 10% between min and max
    minval    = min(min(to_be_plotted));
    maxval    = max(max(to_be_plotted));
    threshold = threshold*(maxval - minval) - minval; % Convert from percentage to absolute value
    is_above_threshold = to_be_plotted > threshold;

    [my, mx] = find(m1 & m2 & m3 & m4 & is_above_threshold);

    % Draw local maxima positions
    hold on;
    plot(xs(mx), ys(my), 'gx');
    hold off;
  end % if

%  % Show position of P2 & P3 hat
%  if (~isempty(P2hat) && ~isempty(P3hat))
%    hold on;
%    X = [period/P2hat; -period/P2hat];
%    Y = [period/P3hat; -period/P3hat];
%    plot(X, Y, 'go', "markersize", 10);
%    hold off;
%  end % if

  % Show filter boundaries
  if (~isempty(filters))
    hold on;
    for n = 1:rows(filters)
      switch filters(n,3)
        case 0 % Horizontal filter
          Xouter = [xmin, xmin;
                    xmax, xmax];
          Xinner = Xouter;
          Youter = [filters(n,1)+filters(n,2), filters(n,1)-filters(n,2);
                    filters(n,1)+filters(n,2), filters(n,1)-filters(n,2)] + shift_DC(2);
          Yinner = [filters(n,1)+filters(n,2)/2, filters(n,1)-filters(n,2)/2;
                    filters(n,1)+filters(n,2)/2, filters(n,1)-filters(n,2)/2] + shift_DC(2);
        case 1 % Vertical filter
          Xouter = [filters(n,1)+filters(n,2), filters(n,1)-filters(n,2);
                    filters(n,1)+filters(n,2), filters(n,1)-filters(n,2)] + shift_DC(1);
          Xinner = [filters(n,1)+filters(n,2)/2, filters(n,1)-filters(n,2)/2;
                    filters(n,1)+filters(n,2)/2, filters(n,1)-filters(n,2)/2] + shift_DC(1);
          Youter = [ymin, ymin;
                    ymax, ymax];
          Yinner = Youter;
      end % switch
      if (filters(n,4) == 0) % filter is not marked for deletion
        plot(Xouter, Youter, 'g--', Xinner, Yinner, 'g');
      else % filter is marked for deletion
        plot(Xouter, Youter, 'r--', Xinner, Yinner, 'r');
      end % if
    end % for
    hold off;
  end % if

  % Draw colorbar
  colorbar('ylabel', clabel_text);
  apply_colormap(fig_no);

  % Set axis labels
  xlabel("2π ν_l P_1");
  ylabel("ν_t P_1");

  % Get/Set axis limits
  if (~set_axes(fig_no))
    get_axes(0,0,fig_no);
  end % if

  if (isempty(get(get(gca, "title"), "string")) && ~isempty(P2hat) && ~isempty(P3hat))
    title_string = sprintf('P2 = P1/%f\nP3 = P1*%f', period/P2hat, P3hat/period);
    title(title_string);
  end % if

end % function

function Xplot_modenv(src, data, fig_no)

  global fig_handles
  global plot_params

  global tdfs

  global period
  global only_visible

  % Switch to/Create 2DFS figure and keep track of the view window
  fig_no = 7;
  first_time = (fig_handles(fig_no) == 0);

  if (first_time)
    create_figure(0,0,fig_no);
  end % if

  figure(fig_handles(fig_no));

  % Assume: tdfs has already been calculated by now

  % Now calculate the inverse 2DFS (=m_{phi,t} in E&S(2002))
  M = ifft2(tdfs);

  % Do singular-value decomposition on m_{phi,t} = m_phi(phi) * m_t(t)
  [U,S,V] = svd(M);
  m_t     = U(:,1) * sqrt(S(1,1));
  m_phi   = V(:,1) * sqrt(S(1,1));

  % Calculate various x-axes
  if (only_visible)
    tmin = max([floor(plot_params(1,1)/period),0]) + 1;
    tmax = tmin + length(m_t) - 1;
  else
    tmin = 1;
    tmax = length(m_t);
  end % if

  pulse_no = [tmin:tmax];
  longitude = [0:(length(m_phi)-1)] / length(m_phi) * 360;

  subplot(2,2,1);
  plot(pulse_no, abs(m_t));
  xlabel("Pulse Number");
  ylabel("Amplitude");

  subplot(2,2,2);
  plot(longitude, abs(m_phi));
  xlabel("Pulse longitude");
  ylabel("Amplitude");

  subplot(2,2,3);
  plot(pulse_no, arg(m_t) * 180/pi, 'x');
  xlabel("Pulse Number");
  ylabel("Phase");

  subplot(2,2,4);
  plot(longitude, arg(m_phi) * 180/pi, 'x');
  xlabel("Pulse longitude");
  ylabel("Phase");

  % Calculate the "signal-to-noise" of the result
  SNR = S(1,1) / (trace(S) - S(1,1)); % <-- Not sure if this is the correct way to do it...

  set(fig_handles(fig_no), "Name", sprintf('Complex Envelope: \"SNR?\" = %f', SNR));

end % function

function Xapply_shift_DC()

  global shift_DC
  global pos_relative
  global tdfs

  % First, a bit of error checking
  if (isempty(shift_DC))
    shift_DC = [0,0];
  end % if

  if (isempty(pos_relative))
    pos_relative = [0,0];
  end % if

  shift_amount = shift_DC + pos_relative;

  % Convert shift_amount to "pixel" units
  pixel_size = [1, 1/rows(tdfs)];
  shift_amount = round(shift_amount ./ pixel_size);

  % And do the shift (as simple as that!)
  tdfs = shift(tdfs, shift_amount(1), 2); % in the x-direction
  tdfs = shift(tdfs, shift_amount(2), 1); % in the y-direction

end % function

function Xapply_tdfs_filters()

  global tdfs
  global filters

  % Get the x and y values
  nxs = columns(tdfs);           nys = rows(tdfs);
  xmin = -floor((nxs-1)/2);      xmax = ceil((nxs-1)/2);
  ymin = -floor((nys-1)/2);      ymax = ceil((nys-1)/2);
  xshift = -xmin;                yshift = -ymin;

  % Shift so that DC is in the centre
  tdfs = shift(tdfs, xshift, 2);
  tdfs = shift(tdfs, yshift);

  xs = [xmin:xmax];          % Units of   v_l*P1/(2*pi)   or   v_l*P1?
  ys = [ymin:ymax] / nys;    % Units v_t * P1

  [X, Y] = meshgrid(xs, ys);

  nfilters = rows(filters);
  for n = 1:nfilters
    switch filters(n,3)
      case 0 % horizontal filter
        transmission = 2*abs(Y - filters(n,1))/filters(n,2) - 1;
      case 1 % vertical filter
        transmission = 2*abs(X - filters(n,1))/filters(n,2) - 1;
    end % switch
    transmission(transmission > 1) = 1;
    transmission(transmission < 0) = 0;
    tdfs .*= transmission;
  end % for

  % Shift back so that DC is at (1,1)
  tdfs = shift(tdfs, -xshift, 2);
  tdfs = shift(tdfs, -yshift);

end % function

function Xreshape_timeseries_into_grid(do_visible = true)

  global timeseries
  global flattened
  global timeseries_grid
  global plot_params

  global period % Assumed to have been set by this point
  global dt

  global only_visible
  global apply_bps

  % Get the appropriate timeseries
  if (apply_bps)
    to_be_plotted = flattened;
  else
    to_be_plotted = timeseries;
  end % if

  % Are we processing just the visible part?
  if (only_visible && do_visible)
    xax = plot_params(1,1:2);
    min_idx = max([floor(xax(1)/dt)+1, 1]);
    max_idx = min([floor(xax(2)/dt)+1, length(timeseries)]); % <-- Check this
    to_be_plotted = to_be_plotted([min_idx:max_idx]);
  end % if

  % Prepare the grid of values
  N         = length(to_be_plotted);
  t         = ([1:N]' - 0.5) * dt; % -0.5 is to avoid some of the more common computer precision errors
  pulse_no  = floor(t/period) + 1;
  phase_bin = floor(mod(t,period)/dt) + 1; % The "floor" here is problematic. It means some of the pulses
                                           % may be shifted by up to dt/2. I see no other way around this.
  accum_subs = [pulse_no, phase_bin];
  timeseries_grid = accumarray(accum_subs, to_be_plotted);

end % function

function Xclear_mask(src, data)

  global profile_mask
  global breakpoint_mask

  profile_mask = [0,0];
  breakpoint_mask(:) = 1;

  flatten();

  set_unsaved_changes(true);

  % Update figures
  get_axes();
  replot();

end % function

function Xselect_mask(src, data)

  global profile_mask

  global fig_handles
  figure(fig_handles(3)); % Assumes figure is already open

  title('Choose low phase for start of mask');
  [x1, y1, button1] = ginput(1);

  title('Choose high phase for end of mask');
  [x2, y2, button2] = ginput(1);

  profile_mask = [x1,x2];
  set_unsaved_changes(true);
  apply_profile_mask();

  % Update figures
  get_axes();
  replot();

end % function

function Xapply_profile_mask()

  global timeseries
  global profile_mask
  global breakpoint_mask

  global period
  global dt

  % Calculate the phase of all the points in the timeseries
  N = length(timeseries);
  t = [0:(N-1)]' * dt;
  phase = mod(t,period)/period;

  % Set the breakpoint_mask to "0" for points "within" the mask
  if (profile_mask(1) <= profile_mask(2))
    mask_idxs = (phase >= profile_mask(1)) & (phase <= profile_mask(2));
  else
    mask_idxs = (phase >= profile_mask(1)) | (phase <= profile_mask(2));
  end % if
  breakpoint_mask = ~mask_idxs;

end % function

function Xcalc_fft()

  global fig_handles

  global timeseries
  global flattened
  global dt
  global period
  global plot_params

  global spectrum_freqs
  global spectrum_vals

  global zeromean
  global zeropad
  global only_visible
  global apply_hamming
  global apply_hanning
  global apply_bps

  % Are we getting the FFT of the original timeseries or the flattened timeseries?
  if (apply_bps)
    to_be_ffted = flattened;
  else
    to_be_ffted = timeseries;
  end % if

  % Are we processing just the visible part?
  if (only_visible)
    xax = plot_params(1,1:2);
    min_idx = max([floor(xax(1)/dt)+1, 1]);
    max_idx = min([floor(xax(2)/dt)+1, length(timeseries)]); % <-- Check this
    to_be_ffted = to_be_ffted([min_idx:max_idx]);
  end % if

  % Apply zeromean
  if (zeromean)
    to_be_ffted = to_be_ffted - mean(to_be_ffted);
  end % if

  % Apply zero-pad
  if (zeropad)
    n = length(to_be_ffted);
    np = period / dt; % Number of bins per pulse
    n_extra_zeros = round(ceil(n/np)*np) - n;
    to_be_ffted = [to_be_ffted; zeros(n_extra_zeros,1)];
  end % if

  % Get the length of the timeseries to be FFT'd
  n  = length(to_be_ffted);

  % Apply Hamming window
  if (apply_hamming)
    to_be_ffted = hamming(n) .* to_be_ffted;
  end % if

  % Apply Hanning window
  if (apply_hanning)
    to_be_ffted = hanning(n) .* to_be_ffted;
  end % if

  % Calculate the values...
  % ...for the FFT abscissa...
  df              = 1/(dt*n);
  spectrum_freqs  = [0:(n-1)] * df;
  % ...and the FFT ordinate (=power)...
  spectrum_vals = fft(to_be_ffted);

end % function

function Xcalc_profile()

  global dt
  global timeseries
  global flattened
  global nprofile_bins
  global profile

  global apply_bps

  global period

  % Use original or flattened, accordingly
  if (apply_bps)
    to_be_folded = flattened;
  else
    to_be_folded = timeseries;
  end % if

  % Calculate the values for the profile abscissa (=phase)
  N = length(to_be_folded);
  t = ([1:N]' - 0.5) * dt; % -0.5 is to avoid some of the more common computer precision errors
  phase = mod(t,period)/period;
  accum_subs = floor(phase * nprofile_bins) + 1;
  profile = accumarray(accum_subs, to_be_folded, [nprofile_bins,1], @mean);

end % function

function Xset_period(newperiod)

  global period

  if (period ~= newperiod)
    set_unsaved_changes(true);
  end % if

  % Change period
  period = newperiod;

  % Update the menu
  update_timeseries_menu();

  % Update the relevant plots
  fig_nos = [3,4,5];
  get_axes(0,0,fig_nos);
  replot(fig_nos);

end % function

function Xset_nbins(newnbins)

  global nprofile_bins

  if (nprofile_bins ~= newnbins)
    set_unsaved_changes(true);
  end % if

  % Change period
  nprofile_bins = newnbins;

  % Update the menu
  update_timeseries_menu();

  % Update the plots
  get_axes(0,0,3);
  replot(3);

end % function

function Xupdate_timeseries_menu()

  global period
  global nprofile_bins
  global zeromean
  global zeropad
  global only_visible
  global apply_hamming
  global apply_hanning
  global apply_bps
  global filename
  global filepath

  global m_fft_window_hamming
  global m_fft_window_hanning
  global m_fft_zeromean
  global m_fft_zeropad
  global m_fft_visible
  global m_bp_apply
  global m_profile_plot
  global m_profile_waterfall
  global m_data_save

  set(m_fft_window_hamming, "checked", on_off(apply_hamming));
  set(m_fft_window_hanning, "checked", on_off(apply_hanning));
  set(m_fft_zeromean,       "checked", on_off(zeromean));
  set(m_fft_zeropad,        "checked", on_off(zeropad));
  set(m_fft_visible,        "checked", on_off(only_visible));
  set(m_bp_apply,           "checked", on_off(apply_bps));
  set(m_profile_plot,       "enable",  on_off(period && nprofile_bins));
  set(m_profile_waterfall,  "enable",  on_off(period));
  set(m_data_save,          "enable",  on_off(filename && filepath));

end % function

function Xget_period_from_user(src, data)

  global period

  cstr = inputdlg({"Enter period in seconds:"}, "Period", 1, {num2str(period, 15)});
  if (~isempty(cstr))
    set_period(str2num(cstr{1}));
  end % if

end % function

function Xget_nbins_from_user(src, data)

  global nprofile_bins

  cstr = inputdlg({"Enter number of bins:"}, "Profile bins", 1, {num2str(nprofile_bins)});

  % Check if user clicked cancel
  if (isempty(cstr))
    return;
  end % if

  % Convert the value to an appropriate number
  try
    input = str2num(cstr{1});
  catch
    errordlg("Unable to convert input to numeric type");
    return
  end % try_catch

  % Check if they put in something other than an integer
  if (mod(input,1) ~= 0)
    input = round(input);
    errordlg(sprintf('Rounding %s to %d', cstr, input));
  end % if

  % Update the value!
  set_nbins(input);

end % function

function Xflatten()

  global dt
  global timeseries
  global flattened
  global breakpoint_mask

  global breakpoints
  global ms  % The slopes of the model detrend lines
  global cs  % The y-intercepts of the model detrend lines

  if (isempty(breakpoints))
    flattened = timeseries;
  else
    flattened = zeros(size(timeseries));
  end % if

  % Book-end the breakpoints with initial and final values
  N = length(timeseries);
  t = [0:(N-1)]' * dt;

  bps = [t(1)-1, breakpoints, t(end)+1];

  % Loop through the breakpoints and detrend each segment
  ms = zeros(length(bps)-1, 1);
  cs = zeros(length(bps)-1, 1);
  for i = 1:(length(bps)-1)
    flatten_idxs  = (t >= bps(i)) & (t < bps(i+1));
    model_idxs    = flatten_idxs & breakpoint_mask;
    fit           = polyfit(t(model_idxs), timeseries(model_idxs), 1);
    ms(i)         = fit(1);
    cs(i)         = fit(2);
    y_orig        = timeseries(flatten_idxs);
    x             = t(flatten_idxs);
    flattened(flatten_idxs) ...
                  = y_orig - (ms(i)*x + cs(i));
  end % for

end % function

function Xset_fft_plot_type(src, data, fig_no, newexp)

  global plot_params

  % Decide if we're changing the y-axis or the c-axis
  switch fig_no
    case 2 % FFT plot
      exp_idx = 14; % yexp column in plot_params
      min_idx = 3;  % ymin    "    "      "
      max_idx = 4;  % ymax    "    "      "
    case {4,6} % HDFS, 2DFS
      exp_idx = 16; % cexp column in plot_params
      min_idx = 7;  % cmin    "    "      "
      max_idx = 8;  % cmax    "    "      "
    otherwise
      return % Currently, cannot change this plot parameter for other types of plots
  end % switch

  oldexp_val = ~isnan(plot_params(fig_no, exp_idx)) + 1; % NaN --> 1;  2 --> 2
  newexp_val = ~isnan(newexp) + 1;

  % Calculate "exponent factor" (used for calculating how to change axis scale)
  exp_fact = newexp_val / oldexp_val;

  % Change the value
  plot_params(fig_no, exp_idx) = newexp;

  % Update figure (with recalculated axes)
  get_axes(0,0,fig_no);
  s = sign(plot_params(fig_no,min_idx:max_idx));
  plot_params(fig_no,min_idx:max_idx) = s .* abs(plot_params(fig_no,min_idx:max_idx) .^ exp_fact);
  replot(fig_no);

end % function

function Xreplot(fig_nos = nan)

  global fig_handles
  global fig_functions

  % The default value of NaN means to do ALL figures
  if (isnan(fig_nos))
    fig_nos = [1:length(fig_handles)];
  end % if

  for fig_no = fig_nos

    % Get the figure handle
    h = fig_handles(fig_no);
    if (~h)
      continue
    end % if

    % Run the plotting function
    fig_functions{fig_no}();

  end % for

end % function

function Xselect_pulsarperiod(src, data)

  global pulsarperiods
  global period

  names = pulsarperiods(:,1);

  [sel, ok] = listdlg("ListString", names, "SelectionMode", "Single", "Name", "Select pulsar");
  if (ok)
    set_period(pulsarperiods{sel,2});
  end % if

end % function

function Xget_axes(src, data, fig_nos = nan)

  global plot_params
  global fig_handles

  % The default value of NaN means to do ALL figures
  if (isnan(fig_nos))
    fig_nos = [1:length(fig_handles)];
  end % if

  for fig_no = fig_nos

    % Get the figure handle
    h = fig_handles(fig_no);
    if (~h)
      continue
    end % if

    % Set up "get axis" functions
    all_lim = {@xlim, @ylim, @zlim, @caxis};

    % Loop through x, y, z, and c axes and apply changes

    for ax_no = 1:4

      % Get the relevant axis function
      thislim = all_lim{ax_no};

      % Get plot_params indices
      maxidx = ax_no * 2;
      minidx = maxidx - 1;

      % Ignore z axis if there IS no z axis
      if ((ax_no == 3) && (length(axis()) ~= 6))
        continue
      end % if

      % Set values to current figure's.
      figure(h);
      ax = thislim();

      % Set unsaved changes flag, if necessary
      if ((plot_params(fig_no, minidx) ~= ax(1)) || ...
          (plot_params(fig_no, maxidx) ~= ax(2)))

        set_unsaved_changes(true);

      end % if

      plot_params(fig_no, minidx) = ax(1);
      plot_params(fig_no, maxidx) = ax(2);

    end % for

  end % for

end % function

function Xsuccess = set_axes(fig_no)

  global plot_params
  global fig_handles

  success = true; % This will get changed to false if there are any NaNs in plot_params

  % Get the figure handle
  h = fig_handles(fig_no);
  if (~h)
    return
  end % if

  % Set up "get axis" functions
  all_lim = {@xlim, @ylim, @zlim, @caxis};

  % Loop through x, y, z, and c axes and apply changes

  for ax_no = 1:4

    % Get the relevant axis function
    thislim = all_lim{ax_no};

    % Get plot_params indices
    maxidx = ax_no * 2;
    minidx = maxidx - 1;

    % Get current values from plot_params
    curr_min = plot_params(fig_no, minidx);
    curr_max = plot_params(fig_no, maxidx);

    % Ignore z axis if there IS no z axis
    if ((ax_no == 3) && (length(axis()) ~= 6))
      continue
    end % if

    % If values are good, apply to current figure.
    % Otherwise, set values to current figure's.
    if (~isnan(curr_min) && ~isnan(curr_max))
      figure(h);
      thislim([curr_min, curr_max]);
    else
      success = false;
    end % if

  end % for

end % function

function Xapply_colormap(fig_no)

  global plot_params
  global fig_handles

  % Switch to the appropriate figure
  h = fig_handles(fig_no);
  if (h)
    figure(h);
  else
    return
  end % if

  % For some reason, the listed colormaps don't change the SIZE
  % of the existing colormap, so I have to do that explicitly:
  colormap(zeros(64,3));
  colormap_name = colormap("list"){plot_params(fig_no, 17)};
  to_be_inverted = plot_params(fig_no, 18);
  cmap = colormap(colormap_name);
  if (to_be_inverted)
    colormap(flipud(cmap));
  end

end % function

function Xset_colormap_type(src, data, fig_no, newtype)

  global plot_params

  curr_type = plot_params(fig_no, 17);

  if (newtype ~= curr_type)
    set_unsaved_changes(true);
  end % if

  plot_params(fig_no, 17) = newtype;

  get_axes();
  replot(fig_no);

end % function

function Xclick_p2p3(src, data, fig_no)

  global fig_handles
  global quantise

  global orig_windowbuttondownfcn_p2p3
  global orig keypressfcn_p2p3

  h = fig_handles(fig_no);

  orig_windowbuttondownfcn_p2p3 = get(h, "windowbuttondownfcn");
  orig_keypressfcn_p2p3         = get(h, "keypressfcn");

  title_stem = ["Choose P2 and P3 (hat) by clicking on the 2DFS\nPress 'q' to turn quantisation"];
  title([title_stem, " ", on_off(~quantise)]);

  set(h, "windowbuttondownfcn", @collect_p2p3_click);
  set(h, "keypressfcn", {@collect_quantise_keypress, title_stem});

end % function

function Xcollect_p2p3_click(src, button)

  global period
  global P2hat
  global P3hat
  global timeseries_grid

  global orig_windowbuttondownfcn_p2p3
  global orig_keypressfcn_p2p3

  global quantise

  if (button == 1) % Left mouse click

    % Get the clicked point
    pos = get_curr_pos(gca, quantise);

    % Calculate P2 and P3 hat (remember to round for pixels)
    P2hat = period/pos(1);
    P3hat = period/pos(2);
    title("");
    set(gcf(), "windowbuttondownfcn", orig_windowbuttondownfcn_p2p3);
    set(gcf(), "keypressfcn", orig_keypressfcn_p2p3);

  end % if

  set_unsaved_changes(true);
  get_axes(0,0,6);
  replot(6);

end % function

function Xcreate_filter(src, data, fig_no, hor_or_vert)
% hor_or_vert should be "horizontal" or "vertical"

  global fig_handles
  global new_filter
  global quantise

  global orig_windowbuttondownfcn_filter
  global orig_windowbuttonmotionfcn_filter
  global orig keypressfcn_filter

  hv = find(strcmp(hor_or_vert, {"horizontal", "vertical"})) - 1; % "h..." --> 0;  "v..." --> 1
  new_filter = [0, 0, hv, 0]; % [coord, width, hor_or_vel, marked_for_deletion]
  h = fig_handles(fig_no);

  orig_windowbuttondownfcn_filter   = get(h, "windowbuttondownfcn");
  orig_windowbuttonmotionfcn_filter = get(h, "windowbuttonmotionfcn");
  orig_keypressfcn_filter           = get(h, "keypressfcn");

  title_stem = ["Click the centre of the ", hor_or_vert, " filter\nPress 'q' to turn quantisation"];
  title([title_stem, " ", on_off(~quantise)]);

  set(h, "windowbuttondownfcn", {@collect_filter_click, 1, hor_or_vert}); % "1" = First click
  set(h, "windowbuttonmotionfcn", []); % <-- YET-TO-IMPLEMENT
  set(h, "keypressfcn", {@collect_quantise_keypress, title_stem});

  get_axes(0,0,fig_no);

end % function

function Xcollect_filter_click(src, button, click_no, hor_or_vert)

  global quantise
  global new_filter
  global filters

  global orig_windowbuttondownfcn_filter
  global orig_windowbuttonmotionfcn_filter
  global orig_keypressfcn_filter

  % Get the last clicked position
  pos = get_curr_pos(gca, quantise);
  idx = 2 - new_filter(3); % "h..." --> 2;  "v..." --> 1

  switch click_no
    case 1
      % Include the clicked value in the new filter
      new_filter(1) = pos(idx);

      % Change figure title
      title_stem = ["Click the edge of the ", hor_or_vert, " filter\nPress 'q' to turn quantisation"];
      title([title_stem, " ", on_off(~quantise)]);

      % Change callback functions on mouse click and key-press
      set(src, "windowbuttondownfcn", {@collect_filter_click, 2, hor_or_vert});
      set(src, "keypressfcn", {@collect_quantise_keypress, title_stem});
    case 2
      % Include the clicked value in the new filter and annex the new filter onto the list of filters
      new_filter(2) = abs(pos(idx) - new_filter(1));
      if (new_filter(2) > 0)
        no_filters = rows(filters);
        filters((no_filters+1),:) = new_filter;
        set_unsaved_changes(true);
      else
        errordlg("Cannot create filter of zero width");
      end % if

      % Reset new_filter to empty (so that nothing extraneous is plotted in plot_tdfs)
      new_filter = [];

      % Set things back to the way they were
      set(src, "windowbuttondownfcn", orig_windowbuttondownfcn_filter);
      set(src, "windowbuttonmotionfcn", orig_windowbuttonmotionfcn_filter);
      set(src, "keypressfcn", orig_keypressfcn_filter);
  end % switch

  replot(6);

end % function

function Xcollect_quantise_keypress(src, evt, newtitle_stem)

  if (evt.Character == 'q')

    global quantise
    quantise = ~quantise;

    title([newtitle_stem, " ", on_off(~quantise)]);

  end

end % function

function Xpos = get_curr_pos(h, quantised)

  point = get(h, "currentpoint");
  point = point(1,1:2);

  % If a figure handle was passed, then the "point" is in "figure"
  % coordinates (i.e. in pixels measured from lower-left corner of
  % plot window, and needs to be converted to axes coordinates
  if (isfigure(h))
    hax       = gca();
    ax_pos    = get(hax, "position");
    fig_pos   = get(h, "position");
    fig_width = fig_pos(3:4);
    ax_lims   = axis()(1:4);
    ax_width  = diff(reshape(ax_lims,2,2));
    ax_orig   = ax_lims(1,3);
    point     = (point./fig_width - ax_pos(1:2)) ./ ax_pos(3:4).*ax_width + ax_orig;
  end % if

  if (quantised)
    global timeseries_grid
    dx = 1;                               % Size of a pixel in the x-direction (in the 2DFS)
    dy = 1/rows(timeseries_grid);         % Size of a pixel in the y-direction (in the 2DFS)
    pos = [round(point(1)/dx) * dx, round(point(2)/dy) * dy];
  else
    pos = [point(1), point(2)];
  end % if

end % function

function Xdelete_filter(src, data, fig_no)

  global fig_handles
  global orig_windowbuttondownfcn_delfilter
  global orig_keypressfcn_delfilter
  global filters

  % "Save" the old
  h = fig_handles(fig_no);
  orig_windowbuttondownfcn_delfilter = get(h, "windowbuttondownfcn");

  set(h, "windowbuttondownfcn", @collect_deletefilter_click);
  set(h, "keypressfcn", @keypress_deletefilter);

  title("Left click = select/deselect\nRight click = cancel\n'd' = delete")

  get_axes(0,0,fig_no);

end % function

function Xcollect_deletefilter_click(src, button)

  global filters
  global shift_DC
  global orig_windowbuttondownfcn_delfilter
  global orig_keypressfcn_delfilter

  switch button
    case 1 % Left click = mark for deletion
      % Get current position of cursor and adjust by shift_DC amount
      pos = get_curr_pos(gca, false) - shift_DC;

      % Select filters containing cursor
      h_distances = (filters(:,3) == 0) .* (abs(pos(2) - filters(:,1)) ./ filters(:,2)); % horizontal filters
      v_distances = (filters(:,3) == 1) .* (abs(pos(1) - filters(:,1)) ./ filters(:,2)); % vertical filters
      distances = h_distances + v_distances;
      [mindist, closest] = min(distances);

      if (mindist <= 1)
        filters(closest,4) = ~filters(closest,4);
        replot(6);
        title("Left click = select/deselect\nRight click = cancel\n'd' = delete")
      end % if

    case 3 % Right click = cancel
      filters(:,4) = 0;
      if (isempty(filters))
        clear filters
      end % if

      set(src, "windowbuttondownfcn", orig_windowbuttondownfcn_delfilter);
      set(src, "keypressfcn", orig_keypressfcn_delfilter);
      title("");
      replot(6);
  end % switch

end % function

function Xkeypress_deletefilter(src, evt)

  global filters
  global orig_windowbuttondownfcn_delfilter
  global orig_keypressfcn_delfilter

  if (evt.Character == 'd')
    to_be_deleted = filters(:,4);
    if (any(to_be_deleted))
      set_unsaved_changes(true);
      filters = filters(~filters(:,4),:);
    end % if

    if (isempty(filters))
      clear filters
    end % if

    set(src, "windowbuttondownfcn", orig_windowbuttondownfcn_delfilter);
    set(src, "keypressfcn", orig_keypressfcn_delfilter);
    replot(6);
  end % if

end % function

function Xclick_shiftDC(src, data, fig_no)

  global fig_handles
  global orig_windowbuttondownfcn_shiftDC

  h = fig_handles(fig_no);
  orig_windowbuttondownfcn_shiftDC = get(h, "windowbuttondownfcn");

  get_axes(fig_no);

  set(h, "windowbuttondownfcn", @buttondown_shiftDC);

  title("Click and drag to shift the origin\nRight click to confirm");

end % function

function Xbuttondown_shiftDC(src, button)

  switch button
    case 1 % Left button = initiate drag
      global pos_start
      global orig_windowbuttonmotionfcn_shiftDC
      global orig_windowbuttonupfcn_shiftDC

      pos_start = get_curr_pos(gcf, true);

      orig_windowbuttonmotionfcn_shiftDC = get(src, "windowbuttonmotionfcn");
      orig_windowbuttonupfcn_shiftDC     = get(src, "windowbuttonupfcn");

      set(src, "windowbuttonmotionfcn", @buttonmotion_shiftDC)
      set(src, "windowbuttonupfcn",     @buttonup_shiftDC);
    case 3 % Right button = confirm and stop
      global orig_windowbuttondownfcn_shiftDC
      set(src, "windowbuttondownfcn", orig_windowbuttondownfcn_shiftDC);
      title("");
  end % switch
  
end % function

function Xbuttonmotion_shiftDC(src, button)

  global pos_start
  global pos_relative % <-- this global variable will get used in plot_tdfs()
                      %     for "live" shifting feedback

  % Get current (quantised) mouse position
  pos_curr = get_curr_pos(gcf, true);

  % Work out how much has been moved since the initial button down
  pos_relative = pos_curr - pos_start;

  % Update plot
  replot(6);
  title("Click and drag to shift the origin\nRight click to confirm");

end % function

function Xbuttonup_shiftDC(src, button)

  global pos_relative
  global shift_DC

  global orig_windowbuttonmotionfcn_shiftDC
  global orig_windowbuttonupfcn_shiftDC

  % Update shift_DC and reset pos_relative
  if (any(pos_relative ~= 0))
    set_unsaved_changes(true);
    shift_DC        += pos_relative;
    pos_relative(:)  = 0;
  end % if

  set(src, "windowbuttonmotionfcn", orig_windowbuttonmotionfcn_shiftDC);
  set(src, "windowbuttonupfcn",     orig_windowbuttonupfcn_shiftDC);

end % function

