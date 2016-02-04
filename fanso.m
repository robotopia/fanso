function fanso()
% -- Function: fanso ()
%     Analyse a timeseries using the Fourier ANalysis Suite for Octave.
%     Requires FLTK.
%
%     Written by Sam McSweeney, 2016, Creative Commons Licence
%     sammy.mcsweeney@gmail.com

  clear -global

  % Initialise global variables relating to plots
  init_plot_variables();

  % Load known pulsar periods
  load_periods();

  % Setup figure handling system
  global fig_handles;
  global fig_functions;
  fig_handles = zeros(10,1); % 10 is arbitrary. Should be more than enough. Increase if needed.
                             % (1) = timeseries plot
                             % (2) = FFT plot
                             % (3) = profile plot
                             % (4) = harmonic resolved fluctuation spectrum plot
                             % (5) = waterfall plot of timeseries modulo period
  fig_functions = {@plot_timeseries, @plot_fft, @plot_profile, @plot_hrfs, @plot_waterfall, @plot_tdfs};

  % Require that we are using the correct graphics toolkit (FLTK)
  gtk = 'fltk';
  if (~any(strcmp(available_graphics_toolkits(), gtk)))
    disp(['This function requires ',gtk,' to be installed.'])
    return;
  else
    graphics_toolkit(gtk);
  end % if

  % Draw the (empty) plot for the first time
  plot_timeseries();

end % function

function init_plot_variables()
% Initialise global variables that need resetting upon program startup
% and whenever a new time series is imported.

  % The "main" global variables relating to the data themselves
  global dt                  % the time between samples in the timeseries (in seconds)
  global timeseries          % the original timeseries
  global flattened           % timeseries flattened by linear detrending between breakpoints
  global timeseries_grid     % timeseries recase into a grid (one pulse per row)
  global breakpoint_mask     % 0 = do not use this value in calculating linear trends; 1 = use this value
  global profile_mask        % A pair of phases that define a region of phases to be ignored in the breakpoint linear fits
  global fft_plot_type       % 0 = amplitudes;  1 = power
  global waterfall_plot_type % 0 = 2D color;  1 = 3D heights
  global period              % The folding period (P1)
  global nprofile_bins       % The number of bins to be used for folding output

  % Variables for the size and position of the figure windows
  global screensize
  global gapsize

  % Variables for the plot settings
  global zeromean       % 0 = Do nothing;               1 = Zero mean before applying FFT
  global zeropad        % 0 = Do nothing;               1 = Zero-pad to a nearly-whole number of periods
  global only_visible   % 0 = FFT of entire timeseries; 1 = FFT of only visible timeseries
  global apply_hamming  % 0 = Do nothing;               1 = Apply Hamming window
  global apply_hanning  % 0 = Do nothing;               1 = Apply Hanning window
  global apply_bps      % 0 = Do nothing;               1 = Apply breakpoints (i.e. "flatten" timeseries)

  % A global variable for the breakpoints
  global breakpoints;

  % A global variable for keeping track of the name of the loaded file
  global filename

  timeseries          = zeros(100,1);
  dt                  = 1;
  flattened           = zeros(size(timeseries));
  breakpoint_mask     = ones(size(timeseries));
  profile_mask        = [0,0];
  fft_plot_type       = 0;
  waterfall_plot_type = 0;
  clear period
  clear nprofile_bins

  screensize = get(0, 'screensize');
  gapsize    = 100;

  zeromean      = 0;
  zeropad       = 0;
  only_visible  = 0;
  apply_hamming = 0;
  apply_hanning = 0;
  apply_bps     = 0;

  breakpoints = [];

  filename = [];

end % function

function save_data(filepathname)

  % Get all relevant global variables in this scope
  global dt
  global timeseries
  global flattened
  global breakpoint_mask
  global profile_mask
  global zeromean
  global zeropad
  global only_visible
  global apply_hamming
  global apply_hanning
  global apply_bps
  global breakpoints
  global nprofile_bins
  global period
  global hrfs_cmin
  global hrfs_cmax
  global waterfall_cmin
  global waterfall_cmax
  global tdfs_cmin
  global tdfs_cmax

  % Save all the info!
  save(filepathname, ...
        "dt", ...
        "timeseries", ...
        "flattened", ...
        "breakpoint_mask", ...
        "profile_mask", ...
        "zeromean", ...
        "zeropad", ...
        "only_visible", ...
        "apply_hamming", ...
        "apply_hanning", ...
        "apply_bps", ...
        "breakpoints", ...
        "nprofile_bins", ...
        "period", ...
        "hrfs_cmin", ...
        "hrfs_cmax", ...
        "waterfall_cmin", ...
        "waterfall_cmax", ...
        "tdfs_cmin", ...
        "tdfs_cmax");

end % function

function save_fan(src, data)

  global filepath
  global filename

  % Assumption: This function is only able to be called
  %             if filepath and filename have proper values.

  save_data([filepath, filename]);

end % function

function saveas_fan(src, data)

  global filepath
  global filename

  % Open up a Save File dialog box
  [savefile, savepath] = uiputfile([filepath, filename]);

  if (~strcmp(savepath, "0")) % then they actually selected a file
    save_data([savepath, savefile]);
    filepath = savepath;
    filename = savename;
  end % if

end % function

function load_fan(src, data)

  % Get all relevant global variables in this scope
  global dt
  global timeseries
  global flattened
  global breakpoint_mask
  global profile_mask
  global zeromean
  global zeropad
  global only_visible
  global apply_hamming
  global apply_hanning
  global apply_bps
  global breakpoints
  global nprofile_bins
  global period
  global hrfs_cmin
  global hrfs_cmax
  global waterfall_cmin
  global waterfall_cmax
  global tdfs_cmin
  global tdfs_cmax

  global filename
  global filepath

  % Open up an Open File dialog box
  [loadfile, loadpath] = uigetfile({"*.fan", "FANSO file"});

  if (~strcmp(loadpath, "0")) % then they actually selected a file
    try
      % Load file contents
      load("-text", [loadpath, loadfile]);

      % Store file name + path in global variables
      filename = loadfile;
      filepath = loadpath;

      % Update menu checks and enables
      update_timeseries_menu();

      % (Re-)plot all
      replot_all();
    catch
      errordlg(["Unable to load file \"", loadfile, "\""]);
    end % try_catch

  end % if

end % function

function import_timeseries(src, data)

  global timeseries

  % Get file from Open File dialog box
  [loadfile, loadpath] = uigetfile();

  if (loadfile ~= 0) % The user has actually chosen a file
    try
      % Load the contents of the selected file into a matrix
      mat = load("-ascii", [loadpath, loadfile]);

      % Reset all the variables
      init_plot_variables();

      % Warn the user if file contains more than one column of numbers
      if (~isvector(mat))
        warndlg("Input file contains multiple columns.\nReading only the first column as timeseries.", ...
                "Multiple columns detected");
      end % if

      % Set the timeseries to the (first) column of values
      timeseries = mat(:,1);

      % Mandate that they supply the sampling rate
      change_dt();

      replot_all();

    catch
      errordlg("This file is in an unreadable format.\nSee the '-ascii' option in Octave's load() function for details", ...
               "Open file error");
    end % try
  end % if

end % function

function change_dt(src, data)

  global dt;

  % Read in the new value via a dialog box
  cstr = inputdlg("Please enter the timestep in seconds:", "Set timestep");

  if (~isempty(cstr))
    dt   = str2num(cstr{1});
    replot_all();
  end % if

end % function

function toggle_hamming(src, data)

  global apply_hamming
  apply_hamming = ~apply_hamming;

  if (apply_hamming)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  replot(2, "none");

end % function

function toggle_hanning(src, data)

  global apply_hanning
  apply_hanning = ~apply_hanning;

  if (apply_hanning)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  replot(2, "none");

end % function

function toggle_zeromean(src, data)

  global zeromean
  zeromean = ~zeromean;

  if (zeromean)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  replot(2, "none");

end % function

function toggle_zeropad(src, data)

  global zeropad
  global period
  global m_fft_zeropad

  zeropad = ~zeropad;

  if (nargin < 2)
    src = m_fft_zeropad;
  end % if

  if (zeropad)
    if (isempty(period))
      get_period_from_user();
    end % if
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  replot(2); % FFT

end % function

function toggle_visible(src, data)

  global only_visible
  only_visible = ~only_visible;

  if (only_visible)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  replot(2, "none");
  replot(4, "none");

end % function

function add_breakpoints(src, data)

  global fig_handles;
  global breakpoints;
  global apply_bps;
  global m_bp_apply;

  if (apply_bps)
    toggle_flatten(m_bp_apply);
  end % if

  do
    figure(fig_handles(1));
    title("Left mouse button = add breakpoint; right = remove breakpoint; 's' = stop\nDON'T CLOSE THIS WINDOW!");
    [x, y, button] = ginput(1);
    switch button
      case 1 % left mouse button
        breakpoints = union(breakpoints, x);
      case 3 % right mouse button
        if (length(breakpoints) > 0)
          [nearest_bpdiff, nearest_idx] = min(abs(breakpoints - x));
          breakpoints = setdiff(breakpoints, breakpoints(nearest_idx));
        end
    end % switch
    replot(1, "none");
  until (button == 115) % 115='s'
  title('');

end % function

function toggle_flatten(src, data)

  global fig_handles; % <-- is this (and other similar lines) needed?

  global apply_bps
  apply_bps = ~apply_bps;

  if (apply_bps)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  replot_all("none");

end % function

function create_figure(src, data, fig_no)

  global fig_handles

  switch fig_no
    case 1 % The timeseries figure
      global screensize
      global gapsize
      global period
      global nprofile_bins

      winsize_x  = screensize(3) - 2*gapsize;
      winsize_y  = floor((screensize(4) - 3*gapsize)/2);
      winpos_x   = gapsize;
      winpos_y   = 2*gapsize + winsize_y;

      fig_handles(fig_no) = figure("Name", "Timeseries", ...
                                   "Position", [winpos_x, winpos_y, winsize_x, winsize_y], ...
                                   "DeleteFcn", {@destroy_figure, fig_no});

      if (isempty(period))
        plot_profile_enable_state = "off";
        plot_waterfall_enable_state = "off";
      else
        plot_profile_enable_state = "on";
        plot_waterfall_enable_state = "on";
      end % if

      if (isempty(nprofile_bins))
        plot_profile_enable_state = "off";
      end % if

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up the menu for the timeseries figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % The Data menu
      global m_data_save
      m_data               = uimenu('label', '&Data');
      m_data_open          = uimenu(m_data, 'label', '&Open', 'callback', @load_fan);
      m_data_save          = uimenu(m_data, 'label', '&Save', 'callback', @save_fan);
      m_data_saveas        = uimenu(m_data, 'label', 'Save &As', 'callback', @saveas_fan);
      m_data_import_ts     = uimenu(m_data, 'label', '&Import timeseries', 'separator', 'on', 'callback', @import_timeseries);
      m_data_change_dt     = uimenu(m_data, 'label', '&Change dt', 'callback', @change_dt);

      % The FFT menu
      global m_fft_window_hamming  m_fft_window_hanning  m_fft_zeromean  m_fft_zeropad  m_fft_visible
      m_fft                = uimenu('label', 'FF&T');
      m_fft_window         = uimenu(m_fft, 'label', '&Windowing function', 'accelerator', 'w');
      m_fft_window_hamming = uimenu(m_fft_window, 'label', 'Ha&mming', 'accelerator', 'm', 'callback', @toggle_hamming);
      m_fft_window_hanning = uimenu(m_fft_window, 'label', 'Ha&nning', 'accelerator', 'n', 'callback', @toggle_hanning);
      m_fft_zeromean       = uimenu(m_fft, 'label', '&Zero-mean', 'accelerator', 'z', 'callback', @toggle_zeromean);
      m_fft_zeropad        = uimenu(m_fft, 'label', '&Zero-pad', 'callback', @toggle_zeropad);
      m_fft_visible        = uimenu(m_fft, 'label', 'Only &visible', 'accelerator', 'v', 'callback', @toggle_visible);
      m_fft_plotfft        = uimenu(m_fft, 'label', 'Plot FFT', 'separator', 'on', 'accelerator', 't', 'callback', @plot_fft);
      m_fft_plothrfs       = uimenu(m_fft, 'label', 'Plot HRFS', 'accelerator', 'h', 'callback', @plot_hrfs);
      m_fft_plottdfs       = uimenu(m_fft, 'label', 'Plot 2DFS', 'accelerator', '2', 'callback', @plot_tdfs);

      % The Breakpoints menu
      global m_bp_apply
      m_bp                 = uimenu('label', '&Breakpoints');
      m_bp_add             = uimenu(m_bp, 'label', 'Add &breakpoints', 'accelerator', 'b', 'callback', @add_breakpoints);
      m_bp_apply           = uimenu(m_bp, 'label', 'A&pply breakpoints', 'separator', 'on', ...
                                          'callback', @toggle_flatten);

      % The Profile menu
      global m_profile_plot  m_profile_waterfall
      m_profile            = uimenu('label', '&Profile');
      m_profile_setperiod  = uimenu(m_profile, 'label', '&Set period', 'callback', @get_period_from_user);
      m_profile_selperiod  = uimenu(m_profile, 'label', 'Se&lect period from pulsar', 'callback', @select_pulsarperiod);
      m_profile_setnbins   = uimenu(m_profile, 'label', '&Set no. profile bins', 'callback', @get_nbins_from_user);
      m_profile_plot       = uimenu(m_profile, 'label', '&Plot profile', 'separator', 'on', ...
                                               'accelerator', 'p', 'callback', @plot_profile);
      m_profile_waterfall  = uimenu(m_profile, 'label', 'Plot &waterfall', 'callback', @plot_waterfall);

      % Set checks and enables correctly
      update_timeseries_menu();

    case 2 % The FFT figure
      global screensize
      global gapsize

      winsize_x  = screensize(3) - 2*gapsize;
      winsize_y  = floor((screensize(4) - 3*gapsize)/2);
      winpos_x   = gapsize;
      winpos_y   = gapsize;

      fig_handles(fig_no) = figure("Name", "Fourier Transform", ...
                                   "Position", [winpos_x, winpos_y, winsize_x, winsize_y], ...
                                   "DeleteFcn", {@destroy_figure, fig_no});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up menu for  FFT figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      m_fftplot           = uimenu('label', 'FF&T');
      m_fftplot_abs       = uimenu(m_fftplot, 'label', 'Absolute value', 'callback', {@set_fft_plot_type, 0});
      m_fftplot_abs       = uimenu(m_fftplot, 'label', 'Power', 'callback', {@set_fft_plot_type, 1});

      m_fftanalyse        = uimenu('label', '&Analyse');
      m_fftanalyse_period = uimenu(m_fftanalyse, 'label', '&Select period (P1)', 'callback', @select_period);

    case 3 % The profile figure
      fig_handles(fig_no) = figure("Name", "Profile", ...
                                   "DeleteFcn", {@destroy_figure, fig_no});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up menu for profile figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      m_mask        = uimenu('label', '&Mask');
      m_mask_clear  = uimenu(m_mask, 'label', '&Clear', 'callback', @clear_mask);
      m_mask_select = uimenu(m_mask, 'label', '&Select', 'callback', @select_mask);

    case 4 % The harmonic resolved fluctuation spectrum figure
      fig_handles(fig_no) = figure("Name", "Harmonic Resolved Fluctuation Spectrum", ...
                                   "DeleteFcn", {@destroy_figure, fig_no});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up menu for HRFS figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      m_hrfs_dynamicrange          = uimenu('label', '&Dynamic range');
      m_hrfs_dynamicrange_change   = uimenu(m_hrfs_dynamicrange, 'label', 'Set dynamic range limits', ...
                                                                 'callback', {@change_dynamic_range, fig_no, 0, 0});
      m_hrfs_dynamicrange_incmax   = uimenu(m_hrfs_dynamicrange, 'label', 'Increase Max by 10% of range', 'accelerator', '+', ...
                                              'separator', 'on', 'callback', {@change_dynamic_range, fig_no, 0,  0.1});
      m_hrfs_dynamicrange_decmax   = uimenu(m_hrfs_dynamicrange, 'label', 'Decrease Max by 10% of range', 'accelerator', '_', ...
                                                                 'callback', {@change_dynamic_range, fig_no, 0, -0.1});
      m_hrfs_dynamicrange_incmin   = uimenu(m_hrfs_dynamicrange, 'label', 'Increase Min by 10% of range', 'accelerator', '=', ...
                                                                 'callback', {@change_dynamic_range, fig_no,  0.1, 0});
      m_hrfs_dynamicrange_decmin   = uimenu(m_hrfs_dynamicrange, 'label', 'Decrease Min by 10% of range', 'accelerator', '-', ...
                                                                 'callback', {@change_dynamic_range, fig_no, -0.1, 0});

    case 5 % The waterfall plot of the timeseries
      fig_handles(fig_no) = figure("Name", "Waterfall plot", ...
                                   "DeleteFcn", {@destroy_figure, fig_no});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Set up menu for waterfall figure %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      m_waterfall           = uimenu('label', '&Waterfall plot');
      m_waterfall_2D        = uimenu(m_waterfall, 'label', '2D', 'accelerator', '2', ...
                                                  'callback', {@set_waterfall_plot_type, 0});
      m_waterfall_3D        = uimenu(m_waterfall, 'label', '3D', 'accelerator', '3', ...
                                                  'callback', {@set_waterfall_plot_type, 1});
      m_waterfall_dynamicrange        = uimenu('label', '&Dynamic range');
      m_waterfall_dynamicrange_change = uimenu(m_waterfall_dynamicrange, ...
                                               'label', 'Set dynamic range limits', ...
                                               'callback', {@change_dynamic_range, fig_no, 0, 0});
      m_waterfall_dynamicrange_incmax = uimenu(m_waterfall_dynamicrange, ...
                                               'label', 'Increase Max by 10% of range', ...
                                               'accelerator', '+', 'separator', 'on', ...
                                               'callback', {@change_dynamic_range, fig_no, 0,  0.1});
      m_waterfall_dynamicrange_decmax = uimenu(m_waterfall_dynamicrange, ...
                                               'label', 'Decrease Max by 10% of range', ...
                                               'accelerator', '_', ...
                                               'callback', {@change_dynamic_range, fig_no, 0, -0.1});
      m_waterfall_dynamicrange_incmin = uimenu(m_waterfall_dynamicrange,
                                               'label', 'Increase Min by 10% of range', ...
                                               'accelerator', '=', ...
                                               'callback', {@change_dynamic_range, fig_no,  0.1, 0});
      m_waterfall_dynamicrange_decmin = uimenu(m_waterfall_dynamicrange, ...
                                               'label', 'Decrease Min by 10% of range', ...
                                               'accelerator', '-', ...
                                               'callback', {@change_dynamic_range, fig_no, -0.1, 0});

    case 6 % The 2DFS plot
      fig_handles(fig_no) = figure("Name", "2D Fluctuation Spectrum", ...
                                   "DeleteFcn", {@destroy_figure, fig_no});

      m_tdfs_dynamicrange        = uimenu('label', '&Dynamic range');
      m_tdfs_dynamicrange_change = uimenu(m_tdfs_dynamicrange, ...
                                          'label', 'Set dynamic range limits', ...
                                          'callback', {@change_dynamic_range, fig_no, 0, 0});
      m_tdfs_dynamicrange_incmax = uimenu(m_tdfs_dynamicrange, ...
                                          'label', 'Increase Max by 10% of range', ...
                                          'accelerator', '+', 'separator', 'on', ...
                                          'callback', {@change_dynamic_range, fig_no, 0,  0.1});
      m_tdfs_dynamicrange_decmax = uimenu(m_tdfs_dynamicrange, ...
                                          'label', 'Decrease Max by 10% of range', ...
                                          'accelerator', '_', ...
                                          'callback', {@change_dynamic_range, fig_no, 0, -0.1});
      m_tdfs_dynamicrange_incmin = uimenu(m_tdfs_dynamicrange,
                                          'label', 'Increase Min by 10% of range', ...
                                          'accelerator', '=', ...
                                          'callback', {@change_dynamic_range, fig_no,  0.1, 0});
      m_tdfs_dynamicrange_decmin = uimenu(m_tdfs_dynamicrange, ...
                                          'label', 'Decrease Min by 10% of range', ...
                                          'accelerator', '-', ...
                                          'callback', {@change_dynamic_range, fig_no, -0.1, 0});

  end % switch

end % function

function set_waterfall_plot_type(src, data, newvalue)

  global waterfall_plot_type

  waterfall_plot_type = newvalue;
  plot_waterfall();

end % function

function change_dynamic_range(src, data, fig_no, minfactor, maxfactor)

  global hrfs_cmin
  global hrfs_cmax

  global waterfall_cmin
  global waterfall_cmax

  global tdfs_cmin
  global tdfs_cmax

  if (all(~[minfactor, maxfactor])) % Putting in zeros for these parameters brings up a dialog box
    switch fig_no
      case 4 % HRFS
        cstrs = inputdlg({"Min:", "Max:"}, "Enter new dynamic range limits", ...
                         1, {num2str(hrfs_cmin), num2str(hrfs_cmax)});
        if (~isempty(cstrs))
          hrfs_cmin = str2num(cstrs{1});
          hrfs_cmax = str2num(cstrs{2});
        end % if
      case 5 % Waterfall plot
        cstrs = inputdlg({"Min:", "Max:"}, "Enter new dynamic range limits", ...
                         1, {num2str(waterfall_cmin), num2str(waterfall_cmax)});
        if (~isempty(cstrs))
          waterfall_cmin = str2num(cstrs{1});
          waterfall_cmax = str2num(cstrs{2});
        end % if
      case 6 % 2DFS
        cstrs = inputdlg({"Min:", "Max:"}, "Enter new dynamic range limits", ...
                         1, {num2str(tdfs_cmin), num2str(tdfs_cmax)});
        if (~isempty(cstrs))
          tdfs_cmin = str2num(cstrs{1});
          tdfs_cmax = str2num(cstrs{2});
        end % if
    end % switch
  else
    switch fig_no
      case 4 % HRFS
        crange = hrfs_cmax - hrfs_cmin;
        hrfs_cmin = hrfs_cmin + crange*minfactor;
        hrfs_cmax = hrfs_cmax + crange*maxfactor;
      case 5 % Waterfall plot
        crange = waterfall_cmax - waterfall_cmin;
        waterfall_cmin = waterfall_cmin + crange*minfactor;
        waterfall_cmax = waterfall_cmax + crange*maxfactor;
      case 6 % 2DFS
        crange = tdfs_cmax - tdfs_cmin;
        tdfs_cmin = tdfs_cmin + crange*minfactor;
        tdfs_cmax = tdfs_cmax + crange*maxfactor;
    end % switch
  end % if

  replot(fig_no, "none");

end % function

function str = on_off(bool)
% Converts 1,0 to "on","off" respectively
  if (bool)
    str = "on";
  else
    str = "off";
  end % if
end % function

function destroy_figure(src, data, fig_no)

  global fig_handles

  if (fig_no == 1)
    % Check if there are unsaved changes
    % ... (yet to implement)

    % Close all (open) figures attached to this instance of FANSO
    for fig_no = 2:length(fig_handles)
      h = fig_handles(fig_no);
      if (h)
        close(h);
      end % if
    end % for
  else
    fig_handles(fig_no) = 0;
  end % if

end % function

function plot_timeseries(src, data)

  global fig_handles

  global filename

  global timeseries
  global flattened
  global dt

  global breakpoints
  global apply_bps

  % Switch to/Create timeseries figure and keep track of the view window
  fig_no = 1;
  first_time = (fig_handles(fig_no) == 0);

  if (first_time)
    create_figure(0,0,fig_no);
  end % if

  % Create figure window title
  if (isempty(filename))
    figure_name = "Timeseries";
  else
    figure_name = ["Timeseries: ", filename];
  end % if

  figure(fig_handles(fig_no), "Name", figure_name);

  % Calculate the values for the timeseries abscissa
  N = length(timeseries);
  t = [0:(N-1)]' * dt;

  % Plot different things depending on whether the breakpoints are
  % to be "applied" or not (i.e. whether the timeseries has been flattened).
  if (apply_bps)
    flatten();
    plot(t, flattened, 'b');
  else
    if(~isempty(breakpoints))
      global ms
      global cs

      % Prepare vertical lines for breakpoint plotting
      bp_xs = [1;1] * breakpoints;
      ymin = min(timeseries);
      ymax = max(timeseries);
      bp_ys = repmat([ymin;ymax],size(breakpoints));

      % Prepare the detrend lines for plotting
      model_xs = [0, breakpoints;
                  breakpoints, t(end)];
      m_mat = repmat(ms',2,1);
      c_mat = repmat(cs',2,1);
      model_ys = model_xs .* m_mat + c_mat;

      % Plot the original timeseries, breakpoints, and detrend lines
      plot(t, timeseries, 'b', ...
           bp_xs, bp_ys, 'r', 'linewidth', 2.0, ...
           model_xs, model_ys, 'g', 'linewidth', 2.0);
    else
      % Plot the timeseries!
      plot(t, timeseries, 'b');
    end % if
  end % if
  %figure(fig_handles(fig_no)); % <-- Possibly redundant, remains to be checked
  xlabel('Time (s)');
  ylabel('Timeseries values');

end % function

function plot_fft(src, data)

  global fig_handles

  global spectrum_freqs
  global spectrum_vals
  global fft_plot_type

  % Switch to/Create FFT figure and keep track of the view window
  fig_no = 2;
  first_time = (fig_handles(fig_no) == 0);

  if (first_time)
    create_figure(0,0,fig_no);
  end % if

  figure(fig_handles(fig_no));

  % Calculate the FFT
  calc_fft();

  % What kind of plot?
  switch fft_plot_type
    case 0 % Plot absolute values
      to_be_plotted = abs(spectrum_vals);
      ylabel_text = 'Amplitude';
    case 1 % Plot power spectrum
      to_be_plotted = abs(spectrum_vals).^2;
      ylabel_text = 'Power';
    otherwise
      error("Unknown FFT plot type");
  end % switch

  % Plot up the FFT
  figure(fig_handles(fig_no));
  plot(spectrum_freqs, to_be_plotted, 'b');
  xlabel('Frequency (Hz)');
  ylabel(ylabel_text);

end % function

function select_period(src, data)

  global fig_handles
  figure(fig_handles(2));

  global period;

  title("Click on a harmonic of the desired frequency");

  [x, y, button] = ginput(1);

  % Get the harmonic number via a dialog box
  cstr = inputdlg({"Harmonic number of selected point:"}, "Harmonic", 1, {"1"});
  if (~isempty(cstr))
    nharm = str2num(cstr{1});
    set_period(nharm/x);
  end % if

  % Reset title
  figure(fig_handles(2));
  title("");

end % function

function plot_profile(src, data)

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

  figure(fig_handles(fig_no));

  % Calculate profile
  calc_profile();

  % Plot it up!
  hold off;
  dphase = 1/nprofile_bins;
  phase = [0:dphase:(1-dphase/2)]; % Ensure that there are the correct number of bins
  plot(phase, profile);
  xlabel('Phase');
  ylabel('Flux (arbitrary units)');
  profile_title = sprintf('Period = %.12f s;    No. of bins = %d', period, nprofile_bins);
  title(profile_title);
  xlim([0,1]);

  % Shade the masked area
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

function plot_hrfs(src, data)

  global fig_handles

  global dt
  global period

  global zeropad

  global spectrum_freqs
  global spectrum_vals

  global hrfs_cmin
  global hrfs_cmax

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

  % Prepare the matrix of values
  to_be_plotted = abs(spectrum_vals);
  N             = length(to_be_plotted);
  nx            = round(N*dt/period);
  ny            = floor(N/(2*nx));
  n             = nx*ny;
  to_be_plotted = to_be_plotted(1:n);
  to_be_plotted = reshape(to_be_plotted,nx,ny)';

  xs = [0:(nx-1)]*spectrum_freqs(2)*period;
  ys = [0:(ny-1)];

  % Plot up the FFT
  figure(fig_handles(fig_no));
  imagesc(xs, ys, to_be_plotted);
  xlabel('Frequency*P');
  ylabel('Harmonic Number');
  axis("xy");
  colormap("gray");
  colorbar('ylabel', 'Amplitude');
  if (~isempty(hrfs_cmin) && ~isempty(hrfs_cmax))
    caxis([hrfs_cmin, hrfs_cmax]);
  end % if

end % function

function plot_waterfall(src, data)

  global fig_handles

  global period
  global timeseries_grid

  global waterfall_plot_type
  global waterfall_cmin
  global waterfall_cmax

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

  reshape_timeseries_into_grid();

  % Get appropriate values for x and y axes
  nxs = columns(timeseries_grid);
  nys =    rows(timeseries_grid);

  xs = [0:(nxs-1)] / nxs;  % Phase
  ys = [1:nys];            % Pulse number

  switch waterfall_plot_type
    case 0 % 2D
      imagesc(xs, ys, timeseries_grid);
      colormap("gray");
      axis("xy");
      colorbar();
      if (~isempty(waterfall_cmin) && ~isempty(waterfall_cmax))
        caxis([waterfall_cmin, waterfall_cmax]);
      end % if
    case 1 % 3D
      [Xs, Ys] = meshgrid(xs, ys);
      waterfall(Xs, Ys, timeseries_grid);
      colormap([0,0,0]); % i.e. all black lines
    otherwise
      error("Unknown waterfall plot type requested");
  end % switch

  xlabel('Phase');
  ylabel('Pulse number');

end % function

function plot_tdfs(src, data)

  global fig_handles

  global timeseries_grid
  global period

  global tdfs_cmin
  global tdfs_cmax

  % Switch to/Create HRFS figure and keep track of the view window
  fig_no = 6;
  first_time = (fig_handles(fig_no) == 0);

  if (first_time)
    create_figure(0,0,fig_no);
  end % if

  figure(fig_handles(fig_no));

  % Calculate the 2DFS
  reshape_timeseries_into_grid();
  tdfs = abs(fft2(timeseries_grid));

  % Chop off the right half (= conjugate of left half for real signal)
  % and rotate vertically so that DC is in centre of grid
  nxs = ceil((columns(tdfs)+1)/2);
  nys = rows(tdfs);
  tdfs = tdfs(:,1:nxs);
  yshift = floor(nys/2);
  tdfs = shift(tdfs, yshift);

  xs = [0:(nxs-1)];                % Units of v_l * P1/(2*pi)
  %xs = [0:(nxs-1)]/period;        % Units of v_l * 2*pi*P1 ??
  ys = [-yshift:(yshift-1)] / nys;  % Units v_t * P1

  imagesc(xs, ys, tdfs);
  colormap("gray");
  axis("xy");
  colorbar();
  if (~isempty(tdfs_cmin) && ~isempty(tdfs_cmax))
    caxis([tdfs_cmin, tdfs_cmax]);
  end % if

  xlabel("2π ν_l P_1");
  ylabel("ν_t P_1");

end % function

function reshape_timeseries_into_grid()

  global timeseries
  global flattened
  global timeseries_grid

  global period % Assumed to have been set by this point
  global dt

  global apply_bps

  % Get the appropriate timeseries
  if (apply_bps)
    to_be_plotted = flattened;
  else
    to_be_plotted = timeseries;
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

function clear_mask(src, data)
  global profile_mask
  global breakpoint_mask

  profile_mask = [0,0];
  breakpoint_mask(:) = 1;

  flatten();

  replot_all();
end % function

function select_mask(src, data)

  global profile_mask

  global fig_handles
  figure(fig_handles(3)); % Assumes figure is already open

  title('Choose low phase for start of mask');
  [x1, y1, button1] = ginput(1);

  title('Choose high phase for end of mask');
  [x2, y2, button2] = ginput(1);

  profile_mask = [x1,x2];
  apply_profile_mask();

  % Update figures
  replot_all();

end % function

function apply_profile_mask()

  global profile_mask
  global breakpoint_mask

  global period
  global dt

  % Calculate the phase of all the points in the timeseries
  N = length(breakpoint_mask);
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

function calc_fft()

  global fig_handles

  global timeseries
  global flattened
  global dt
  global period

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
    figure(fig_handles(1));
    ax1 = axis();
    min_idx = max([floor(ax1(1)/dt)+1, 1]);
    max_idx = min([floor(ax1(2)/dt)+1, length(timeseries)]); % <-- Check this
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

function calc_profile()

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

function set_period(newperiod)

  % Change period
  global period
  period = newperiod;

  % Update the menu
  update_timeseries_menu();

  % Update the plots
  replot(3);
  replot(4);
  replot(5);

end % function

function set_nbins(newnbins)

  % Change period
  global nprofile_bins
  nprofile_bins = newnbins;

  % Update the menu
  update_timeseries_menu();

  % Update the plots
  replot(3);

end % function

function update_timeseries_menu()

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

function get_period_from_user(src, data)

  global period

  cstr = inputdlg({"Enter period in seconds:"}, "Period", 1, {num2str(period)});
  if (~isempty(cstr))
    set_period(str2num(cstr{1}));
  end % if

end % function

function get_nbins_from_user(src, data)

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

function flatten()

  global dt
  global timeseries
  global flattened
  global breakpoint_mask

  global breakpoints
  global ms  % The slopes of the model detrend lines
  global cs  % The y-intercepts of the model detrend lines

  if (isempty(breakpoints))
    flattened = timeseries;
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

function set_fft_plot_type(src, data, newvalue)

  global fft_plot_type
  fft_plot_type = newvalue;

  replot(2, "y");

end % function

function replot(fig_no, rescale = "xy")
% function: replot(fig_no, rescale = "none")
%
% rescale can be "none", "x", "y", or "xy"
%

  global fig_handles
  global fig_functions

  % Get the handle for the figure
  h = fig_handles(fig_no);

  % If the figure already exists and has an associated plot function
  if ((fig_no <= length(fig_functions)) && (h ~= 0))

    % Save the old axis info
    ax = axis();

    % Run the plotting function
    fig_functions{fig_no}();

    % Update the axis
    figure(h);
    switch rescale
      case "none"
        axis(ax);
      case "x"
        ylim([ax(3), ax(4)]);
      case "y"
        xlim([ax(1), ax(2)]);
      % case "xy" happens naturally
    end % switch

  end % if

end % function

function replot_all(rescale = "xy")

  global fig_handles

  % Loop through the figure numbers
  for fig_no = 1:length(fig_handles)
    replot(fig_no, rescale);
  end % for

end % function

function load_periods()

  global pulsarperiods

  f = fopen('periods.dat');

  % First, count the lines
  nlines = fskipl(f, Inf);
  frewind(f);

  % Now read them in
  pulsarperiods = cell(nlines, 2);
  for n = 1:nlines
    [pulsarperiods{n,1}, pulsarperiods{n,2}] = fscanf(f, '%s %f', "C");
  end % for

end % function

function select_pulsarperiod(src, data)

  global pulsarperiods
  global period

  names = pulsarperiods(:,1);

  [sel, ok] = listdlg("ListString", names, "SelectionMode", "Single", "Name", "Select pulsar");
  if (ok)
    set_period(pulsarperiods{sel,2});
  end % if

end % function
