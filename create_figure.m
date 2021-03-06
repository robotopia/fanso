function create_figure(plot_name)

  global data;
  global figures;

  % Create new figure if it doesn't exist. If it does exist, do nothing.
  if (isempty(figures.(plot_name).fig_handle))

    % Get a new figure handle
    f = figure("name", plot_name);
    figures.(plot_name).fig_handle = f;

    % Set to default position
    if (isfield(figures.(plot_name), "defaultpos"))
      set(f, "position", figures.(plot_name).defaultpos);
    end % if

    % Set callback functions for mouse and keyboard events
    if (~isempty(data))
      set(figures.(plot_name).fig_handle, "keypressfcn", @keypressfcn);
      set(figures.(plot_name).fig_handle, "windowbuttondownfcn", {@panzoom_down, plot_name});
    end % if

    % Set callback function for when figure is closed
    set(f, "closerequestfcn", {@close_figure, plot_name});

    % Create a data menu for the timeseries plot
    if (strcmp(plot_name, "timeseries"))
      % Are there data? If not, disable some menu items
      if (isempty(data))
        enable_state = "off";
      else
        enable_state = "on";
      end % if

      global m_data_save m_data_saveas m_data_export
      m_data         = uimenu("label", "&Data");
      m_data_new     = uimenu(m_data, "label", "&New",        "callback", @new_fan);
      m_data_open    = uimenu(m_data, "label", "&Open..",    "callback", @load_fan);
      m_data_save    = uimenu(m_data, "label", "&Save",       "callback", @save_fan,   "enable", enable_state);
      m_data_saveas  = uimenu(m_data, "label", "Save &As..", "callback", @saveas_fan, "enable", enable_state);
      m_data_import  = uimenu(m_data, "label", "&Import timeseries..", "callback", @import_timeseries, "separator", "on");
      m_data_importp = uimenu(m_data, "label", "Import &PRESTO .dat..", "callback", @import_presto_dat);
      m_data_export  = uimenu(m_data, "label", "&Export timeseries..", "callback", @export_timeseries, "enable", enable_state);
    end % if

  end % if

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback functions for the menu items %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function new_fan()

  global figures;
  global data;
  global filename;

  % Ask about unsaved changes
  if (~offer_to_save())
    return
  end % if

  filename = [];
  data = [];

  f = figures.timeseries.fig_handle;
  a = figures.timeseries.ax_handle;

  delete(a);
  set(f, "name", "timeseries");

  figures.timeseries.ax_handle = [];

  % Set enables on menu items
  global m_data_save m_data_saveas m_data_export;
  set([m_data_save, m_data_saveas, m_data_export], "enable", "off");

  % Set callback function for when a key is pressed
  set(figures.timeseries.fig_handle, "keypressfcn", []);
  set(figures.timeseries.fig_handle, "windowbuttondownfcn", []);
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function load_fan()

  global figures;

  % Ask about unsaved changes
  if (~offer_to_save())
    return
  end % if

  % Open up an Open File dialog box
  [loadfile, loadpath] = uigetfile({"*.fan", "FANSO file"});

  if (~strcmp(loadpath, "0")) % then they actually selected a file
    load_data(loadpath, loadfile);

    % Reset unsaved changes flag
    set_unsaved_changes(false);

    % Redraw all open figures
    plot_names = fieldnames(figures);
    for n = 1:length(plot_names)
      plot_name = plot_names{n};
      figures.(plot_name).drawfcn();
    end % if

    % Change which menu items are enabled
    global m_data_save m_data_saveas m_data_export;
    set([m_data_save, m_data_saveas, m_data_export], "enable", "on");

    % Set callback function for when a key is pressed
    % (in case it was unset by new_fan())
    set(figures.timeseries.fig_handle, "keypressfcn", @keypressfcn);
    set(figures.timeseries.fig_handle, "windowbuttondownfcn", {@panzoom_down, "timeseries"});
  end % if

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_fan()

  global filepath;
  global filename;

  % Assumption: This function is only able to be called
  %             if filepath and filename have proper values.

  save_data([filepath, filename]);

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cont = saveas_fan()

  global filepath;
  global filename;
  cont = true;

  % Open up a Save File dialog box
  if (filepath && filename)
    [savefile, savepath] = uiputfile([filepath, filename]);
  elseif (filename)
    [savefile, savepath] = uiputfile(filename);
  else
    [savefile, savepath] = uiputfile({"*.fan", "FANSO file"});
  end % if

  if (~strcmp(savepath, "0")) % then they actually selected a file
    save_data([savepath, savefile]);
    filepath = savepath;
    filename = savefile;

    % Change which menu items are enabled
    global m_data_save
    set(m_data_save, "enable", "on");
  else
    cont = false;
  end % if

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function import_presto_dat()

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  global filename;
  global filepath;

  % Ask about unsaved changes
  if (~offer_to_save())
    return
  end % if

  % Get file from Open File dialog box
  if (isempty(filepath))
    filepath = pwd();
  end % if

  [loadfile, loadpath] = uigetfile([filepath, "/*"]);

  if (loadfile ~= 0) % The user has actually chosen a file
    try
      % Load the contents of the selected file into a vector
      f = fopen([loadpath, loadfile]);
      dat = fread(f, Inf, "float32");
      fclose(f);
    catch
      errordlg("This file is in an unreadable format.\nSee the '-ascii' option in Octave's load() function for details", ...
               "Open file error");
      return
    end % try

    % Delete all figures other than timeseries
    plot_names = fieldnames(figures);
    for n = 1:length(plot_names)

      plot_name = plot_names{n};

      if (strcmp(plot_name, "timeseries"))
        continue;
      end % if

      close(figures.(plot_name).fig_handle);

    end % for

    % Set all the data variables
    data = struct();
    data.timeseries    = dat(:,1);
    data.samplingrate  = 1;
    data.timeunits     = "sec";
    data.frequnits     = "Hz";

    % Calculate the time series abcissa, etc.
    analysed = struct();
    analysed.N  = length(data.timeseries);             % The length of the timeseries
    analysed.dt = 1/data.samplingrate;                 % The time between adjacent samples
    analysed.t  = [0:(analysed.N-1)]' * analysed.dt;   % The time axis

    % Reset all the other "temporary" variables
    plots = struct();

    for n = 1:length(plot_names)
      if (figures.(plot_names{n}).dims == 2)
        plots.(plot_names{n}).cmap        = 8;              % Colormap idx (8 = greyscale)
        plots.(plot_names{n}).cinv        = false;          % Invert colormap?
        plots.(plot_names{n}).cax         = [];             % axis limits for the colormap
      end % if
      if (figures.(plot_names{n}).isfft)
        plots.(plot_names{n}).islog       = 0;              % Display values in linear(=0) | log(=1) scale
        plots.(plot_names{n}).ispower     = 0;              % Display values as amplitudes(=0) | powers(=1)
      end % if
      plots.(plot_names{n}).axis          = [];             % xlim, ylim, etc. For saving to file
    end % for

    % "analysis" structure
    analysis = struct();
    analysis.profile_mask         = []; % A pair of phases that define a region of phases to be ignored in the breakpoint linear fits
    analysis.period               = []; % The folding period
    analysis.P2hat                = []; % The measured longitudinal "time" between subpulses, P2
    analysis.P3hat                = []; % The measured "time" between subpulses at the same phase, P3
    analysis.nprofile_bins        = []; % The number of bins to be used for folding output
    analysis.zeromean             = 0;  % 0 = Do nothing;                1 = Zero mean before applying FFT
    analysis.zeropad              = 0;  % 0 = Do nothing;                1 = Zero-pad to a nearly-whole number of periods
    analysis.only_visible         = 0;  % 0 = Analyse entire timeseries; 1 = Analyse only visible timeseries
    analysis.only_visible_stack   = 0;  % 0 = Analyse entire pulsestack; 1 = Analyse only visible pulsestack
    analysis.apply_hamming        = 0;  % 0 = Do nothing;                1 = Apply Hamming window
    analysis.apply_hanning        = 0;  % 0 = Do nothing;                1 = Apply Hanning window
    analysis.apply_bps            = 0;  % 0 = Do nothing;                1 = Apply breakpoints (i.e. "flatten" timeseries)
    analysis.breakpoints          = []; % The x-coordinates of the breakpoints
    analysis.filters              = []; % Horizontal and vertical filters (E&S 2002 - style) used on the 2DFS
    analysis.shift_DC             = []; % Horizontal and vertical displacement of the origin in the 2DFS

    % Draw the timeseries plot
    drawfcn = figures.timeseries.drawfcn;
    drawfcn();

    % Change which menu items are enabled
    global m_data_save m_data_saveas m_data_export;
    set([m_data_saveas, m_data_export], "enable", "on");
    set(m_data_save, "enable", "off");

    % Set callback function for when a key is pressed
    % (in case it was unset by new_fan())
    set(figures.timeseries.fig_handle, "keypressfcn", @keypressfcn);
    set(figures.timeseries.fig_handle, "windowbuttondownfcn", {@panzoom_down, "timeseries"});

    % Clear the (saved) filename and path variables
    filename = [];
    filepath = loadpath;

    % Are there changes? Yes!
    set_unsaved_changes(true);

  end % if


end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function import_timeseries()

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  global filename;
  global filepath;

  % Ask about unsaved changes
  if (~offer_to_save())
    return
  end % if

  % Get file from Open File dialog box
  if (isempty(filepath))
    filepath = pwd();
  end % if

  [loadfile, loadpath] = uigetfile([filepath, "/*"]);

  if (loadfile ~= 0) % The user has actually chosen a file
    try
      % Load the contents of the selected file into a matrix
      mat = load("-ascii", [loadpath, loadfile]);
    catch
      errordlg("This file is in an unreadable format.\nSee the '-ascii' option in Octave's load() function for details", ...
               "Open file error");
      return
    end % try

    % Warn the user if file contains more than one column of numbers
    if (~isvector(mat))
      warndlg("Input file contains multiple columns.\nReading only the first column as timeseries.", ...
              "Multiple columns detected");
    end % if

    % Delete all figures other than timeseries
    plot_names = fieldnames(figures);
    for n = 1:length(plot_names)

      plot_name = plot_names{n};

      if (strcmp(plot_name, "timeseries"))
        continue;
      end % if

      close(figures.(plot_name).fig_handle);

    end % for

    % Set all the data variables
    data = struct();
    data.timeseries    = mat(:,1);
    data.samplingrate  = 1;
    data.timeunits     = "sec";
    data.frequnits     = "Hz";

    % Calculate the time series abcissa, etc.
    analysed = struct();
    analysed.N  = length(data.timeseries);             % The length of the timeseries
    analysed.dt = 1/data.samplingrate;                 % The time between adjacent samples
    analysed.t  = [0:(analysed.N-1)]' * analysed.dt;   % The time axis

    % Reset all the other "temporary" variables
    plots = struct();

    for n = 1:length(plot_names)
      if (figures.(plot_names{n}).dims == 2)
        plots.(plot_names{n}).cmap        = 8;              % Colormap idx (8 = greyscale)
        plots.(plot_names{n}).cinv        = false;          % Invert colormap?
        plots.(plot_names{n}).cax         = [];             % axis limits for the colormap
      end % if
      if (figures.(plot_names{n}).isfft)
        plots.(plot_names{n}).islog       = 0;              % Display values in linear(=0) | log(=1) scale
        plots.(plot_names{n}).ispower     = 0;              % Display values as amplitudes(=0) | powers(=1)
      end % if
      plots.(plot_names{n}).axis          = [];             % xlim, ylim, etc. For saving to file
    end % for

    % "analysis" structure
    analysis = struct();
    analysis.profile_mask         = []; % A pair of phases that define a region of phases to be ignored in the breakpoint linear fits
    analysis.period               = []; % The folding period
    analysis.P2hat                = []; % The measured longitudinal "time" between subpulses, P2
    analysis.P3hat                = []; % The measured "time" between subpulses at the same phase, P3
    analysis.nprofile_bins        = []; % The number of bins to be used for folding output
    analysis.zeromean             = 0;  % 0 = Do nothing;                1 = Zero mean before applying FFT
    analysis.zeropad              = 0;  % 0 = Do nothing;                1 = Zero-pad to a nearly-whole number of periods
    analysis.only_visible         = 0;  % 0 = Analyse entire timeseries; 1 = Analyse only visible timeseries
    analysis.only_visible_stack   = 0;  % 0 = Analyse entire pulsestack; 1 = Analyse only visible pulsestack
    analysis.apply_hamming        = 0;  % 0 = Do nothing;                1 = Apply Hamming window
    analysis.apply_hanning        = 0;  % 0 = Do nothing;                1 = Apply Hanning window
    analysis.apply_bps            = 0;  % 0 = Do nothing;                1 = Apply breakpoints (i.e. "flatten" timeseries)
    analysis.breakpoints          = []; % The x-coordinates of the breakpoints
    analysis.filters              = []; % Horizontal and vertical filters (E&S 2002 - style) used on the 2DFS
    analysis.shift_DC             = []; % Horizontal and vertical displacement of the origin in the 2DFS

    % Draw the timeseries plot
    drawfcn = figures.timeseries.drawfcn;
    drawfcn();

    % Change which menu items are enabled
    global m_data_save m_data_saveas m_data_export;
    set([m_data_saveas, m_data_export], "enable", "on");
    set(m_data_save, "enable", "off");

    % Set callback function for when a key is pressed
    % (in case it was unset by new_fan())
    set(figures.timeseries.fig_handle, "keypressfcn", @keypressfcn);
    set(figures.timeseries.fig_handle, "windowbuttondownfcn", {@panzoom_down, "timeseries"});

    % Clear the (saved) filename and path variables
    filename = [];
    filepath = loadpath;

    % Are there changes? Yes!
    set_unsaved_changes(true);

  end % if


end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function export_timeseries()
  % YET TO BE IMPLEMENTED
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_data(file)

  global data;
  global plots;
  global analysis;

  save("-binary", file, "data", "plots", "analysis");
  set_unsaved_changes(false);

end % function

function cont = offer_to_save()
% Returns true if either nothing needs saving
% or user chooses either Yes or No from the
% dialog box. Returns false if they choose Cancel.

  global filename;
  global unsaved_changes;

  cont = true;

  if (unsaved_changes)
    dlg_answer = questdlg("Would you like to save your changes?", "Save changes?");

    switch dlg_answer

      case "Yes"
        if (filename)
          save_fan();
        else
          if (~saveas_fan())
            cont = false;
          end % if
        end % if

      case "Cancel"
        cont = false;

    end % switch
  end % if

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function close_figure(src, data, plot_name)

  global figures;

  % The timeseries figure gets special treatment
  if (strcmp(plot_name, "timeseries"))

    % Check if there are unsaved changes
    if (~offer_to_save())
      return
    end % if

    % Close all (open) figures attached to this instance of FANSO
    plot_names = fieldnames(figures);
    for n = 1:length(plot_names)
      plot_name = plot_names{n};
      delete(figures.(plot_name).fig_handle);
    end % if

    % Get rid of all global variables
    clear -global

  else

    % Close just this figure and set the handle to empty
    delete(src);
    figures.(plot_name).fig_handle = [];
    figures.(plot_name).ax_handle  = [];

  end % if

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keypressfcn(src, evt)

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  % Get the plot_name of the figure window in which the key was pressed
  plot_name  = fighandle2plotname(src);
  plot_names = fieldnames(figures);

  if (isempty(plot_name))
    error("Unknown callback figure handle in keypressfcn");
  end % if

  if (any(strcmp(evt.Key, {"f1", "f2", "f3", "f4", "f5", "f6"})))
    plot_no         = str2num(evt.Key(2));
    this_plot_name  = plot_names{plot_no};
    if (isempty(figures.(this_plot_name).fig_handle))
      create_figure(this_plot_name);
      figures.(this_plot_name).drawfcn();
    end % if
    return
  end % switch

  if (strcmp(plot_name, "modenv"))
    if (isfield(analysed, "Sn"))
      switch evt.Key
        case "rightarrow"
          analysed.Sn = min([analysed.Sn+1, analysed.SN]);
          figures.(plot_name).drawfcn();
          return
        case "leftarrow"
          analysed.Sn = max([analysed.Sn-1, 1]);
          figures.(plot_name).drawfcn();
          return
      end % switch
    end  % if
  end % if

  switch evt.Character
    case 'h'
      %%%%%%%%%%%%%%%%
      % Display help %
      %%%%%%%%%%%%%%%%

      msgbox({"F1 = plot timeseries",
              "F2 = plot FFT",
              "F3 = plot profile",
              "F4 = plot pulse stack",
              "F5 = plot 2DFS",
              "F6 = plot modulation envelopes",
              "0 = toggle zero-padding",
              "2 = toggle power / amplitude",
              "+ = increase dynamic range max",
              "= = decrease dynamic range max",
              "_ = increase dynamic range min",
              "- = decrease dynamic range min",
              "b = change number of profile bins",
              "B = toggle breakpoint edit mode",
              "c = change colormap",
              "d = shift DC (in 2DFS)",
              "e = add filter (to 2DFS)",
              "E = delete filter (from 2DFS)",
              "f = flatten timeseries",
              "h = display this help",
              "i = invert colormap",
              "k = select profile mask",
              "K = clear profile mask",
              "l = toggle logarithmic plot for current figure",
              "m = toggle Hamming window",
              "n = toggle Hanning window",
              "o = change folding period by clicking on FFT",
              "p = change folding period manually",
              "P = change the folding period by providing an ephemeris file",
              "s = change sampling rate",
              "v = toggle only analyse visible part of timeseries",
              "V = toggle only analyse visible part of pulse stack",
              "z = toggle zero-meaning"},
              "Keyboard shortcuts");
    case '0'
      %%%%%%%%%%%%%%%%%%
      % Toggle zeropad %
      %%%%%%%%%%%%%%%%%%
      if (isempty(analysis.period))
        errordlg("Cannot zeropad without first setting period (P1)");
      else
        toggle_analysis_value("zeropad");
        figures.fft.drawfcn();
      end % if

    case '2'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Toggle power / amplitude %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (figures.(plot_name).isfft)
        plots.(plot_name).ispower = ~plots.(plot_name).ispower;
        set_unsaved_changes(true);
        figures.(plot_name).drawfcn();
      end % if

    case '+'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Increase dynamic range max %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (figures.(plot_name).dims == 2) % 2D plot
        cmax = plots.(plot_name).cax(2);
        cmin = plots.(plot_name).cax(1);
        plots.(plot_name).cax(2) += 0.1*(cmax - cmin);
        set_unsaved_changes(true);
        caxis(figures.(plot_name).ax_handle, plots.(plot_name).cax);
      end % if

    case '='
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Decrease dynamic range max %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (figures.(plot_name).dims == 2) % 2D plot
        cmax = plots.(plot_name).cax(2);
        cmin = plots.(plot_name).cax(1);
        plots.(plot_name).cax(2) -= 0.1*(cmax - cmin);
        set_unsaved_changes(true);
        caxis(figures.(plot_name).ax_handle, plots.(plot_name).cax);
      end % if

    case '_'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Increase dynamic range min %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (figures.(plot_name).dims == 2) % 2D plot
        cmax = plots.(plot_name).cax(2);
        cmin = plots.(plot_name).cax(1);
        plots.(plot_name).cax(1) += 0.1*(cmax - cmin);
        set_unsaved_changes(true);
        caxis(figures.(plot_name).ax_handle, plots.(plot_name).cax);
      end % if

    case '-'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Decrease dynamic range min %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (figures.(plot_name).dims == 2) % 2D plot
        cmax = plots.(plot_name).cax(2);
        cmin = plots.(plot_name).cax(1);
        plots.(plot_name).cax(1) -= 0.1*(cmax - cmin);
        set_unsaved_changes(true);
        caxis(figures.(plot_name).ax_handle, plots.(plot_name).cax);
      end % if

    case 'a'
      %%%%%%%%%%%%%%%%%%%%
      % Autoscale figure %
      %%%%%%%%%%%%%%%%%%%%
      if (~isequal(plots.(plot_name).axis, plots.(plot_name).autoscale))
        plots.(plot_name).axis = plots.(plot_name).autoscale;
        set_unsaved_changes(true);

        % Don't replot the figure, just set the axes
        for n = 1:figures.(plot_name).nplots
          axis(figures.(plot_name).ax_handle(n), plots.(plot_name).axis(n,:));
        end % for

        % If they've moved the timeseries plot and only_visible is turned on,
        % then replot everything.
        if (strcmp(plot_name, "timeseries") && analysis.only_visible)
          for n = 1:length(plot_names)
            figures.(plot_names{n}).drawfcn();
          end % for
        end % if

        % If they've moved the pulsestack plot and only_visible_stack is turned on,
        % then replot 2DFS and modulation envelope plots
        if (strcmp(plot_name, "pulsestack") && analysis.only_visible_stack)
          figures.tdfs.drawfcn();
          figures.modenv.drawfcn();
        end % if

      end % if

    case 'b'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Change number of profile bins manually %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      cstr = inputdlg({"Enter number of bins to use in profile:"}, "Profile bins", 1, {num2str(analysis.nprofile_bins)});
      if (~isempty(cstr)) % = OK was pushed
        % Convert the input value to a number
        try
          input = str2num(cstr{1});
        catch
          errordlg("Unable to convert input to numeric type");
          return
        end % try_catch

        % Check if they put in something other than an integer
        if (mod(input,1) ~= 0)
          input = round(input);
          errordlg(sprintf("Rounding %s to %d", cstr, input));
        end % if

        % Update the value!
        set_analysis_value("nprofile_bins", input);

        % Update just the profile figure
        figures.profile.drawfcn();
      end % if

    case 'B'
      %%%%%%%%%%%%%%%%%%%
      % Add breakpoints %
      %%%%%%%%%%%%%%%%%%%
      if (~isfield(analysed, "bp_editmode"))
        analysed.bp_editmode = false;
      end % if

      if (analysed.bp_editmode) % i.e. it is in edit mode when b was pressed
        
        analysis.apply_bps = analysed.old_apply_bps; % i.e. return apply_bps to its former value
        figures.timeseries.drawfcn();

        % If there have been changes, update all figures (they will all be affected)
        if (~isequal(analysed.old_breakpoints, analysis.breakpoints))
          set_unsaved_changes(true);
          plot_names = fieldnames(figures);
          for n = 1:length(plot_names)
            if (~strcmp(plot_names{n}, "timeseries"))
              figures.(plot_names{n}).drawfcn();
            end % if
          end % for
        end % if

        % Clear the title
        analysed = rmfield(analysed, "timeseries_title");
        set_title("timeseries");

        % Reset callback functions
        set(figures.timeseries.fig_handle, "windowbuttondownfcn", []);

      else % i.e. it is NOT in edit mode when b was pressed

        analysed.old_apply_bps    = analysis.apply_bps;
        analysed.old_breakpoints  = analysis.breakpoints;
        analysed.timeseries_title = "Left mouse button = add breakpoint; right = remove breakpoint";

        if (analysis.apply_bps)
          analysis.apply_bps = 0; % i.e. turn it OFF (temporarily)
          figures.timeseries.drawfcn();
        end % if

        set_title("timeseries");

        set(figures.timeseries.fig_handle, "windowbuttondownfcn", @collect_breakpoint_clicks);

      end % if

      % Flip the mode
      analysed.bp_editmode = ~analysed.bp_editmode;

    case 'c'
      %%%%%%%%%%%%%%%%%%%
      % Change colormap %
      %%%%%%%%%%%%%%%%%%%
      [sel, ok] = listdlg("liststring",    colormap("list"), ...
                          "selectionmode", "single", ...
                          "initialvalue",  plots.(plot_name).cmap);

      if (ok)
        plots.(plot_name).cmap = sel;
        figures.(plot_name).drawfcn();
      end % if

    case 'd'
      if (strcmp(plot_name, "tdfs")) % Pressing 'd' only works on 2DFS plot
        if (isfield(analysed, "filter")) % This only happens when in "create filter" mode
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Toggle direction of filter when adding filters %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          switch analysed.filter.direction
            case "horizontal"
              analysed.filter.direction = "vertical";
            case "vertical"
              analysed.filter.direction = "horizontal";
          end % switch
          analysed.tdfs_title = ["Click on the centre of the new (", analysed.filter.direction, ...
                                 ") filter\nPress 'd' to toggle direction filter\n", ...
                                 "Press 'q' to toggle quantisation (currently ", ...
                                 analysed.filter.quantised, ")\n", ...
                                 "Right click to cancel"];
          set_title(plot_name);
        else % When not in any special mode
          %%%%%%%%%%%%
          % Shift DC %
          %%%%%%%%%%%%
          analysed.tdfs_title = "Click on shifted origin\nRight click to cancel";
          set_title(plot_name);
          set(figures.(plot_name).fig_handle, "windowbuttondownfcn", @shift_DC);
        end % if
      end % if

    case 'e'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Create filter (in 2DFS) %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (strcmp(plot_name, "tdfs"))
        analysed.filter.direction = "horizontal"; % h = horizontal
        analysed.filter.quantised = "on";
        analysed.filter.click_no = 1;
        analysed.tdfs_title = ["Click on the centre of the new (", analysed.filter.direction, ...
                               ") filter\nPress 'd' to toggle direction filter\n", ...
                               "Press 'q' to toggle quantisation (currently ", ...
                               analysed.filter.quantised, ")\n", ...
                               "Right click to cancel"];
        figures.(plot_name).drawfcn();
        figures.modenv.drawfcn();
        set(figures.(plot_name).fig_handle, "windowbuttondownfcn", @add_filter);
      end % if

    case 'E'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Delete filter (from 2DFS) %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (strcmp(plot_name, "tdfs") && ~isempty(analysis.filters))
        analysed.filter.editmode = true;
        analysed.tdfs_title = "Click on filter to remove\nRight click to cancel";
        figures.(plot_name).drawfcn();
        figures.modenv.drawfcn();
        set(figures.(plot_name).fig_handle, "windowbuttondownfcn", @remove_filter);
      end % if

    case 'f'
      %%%%%%%%%%%%%%%%%%
      % Toggle flatten %
      %%%%%%%%%%%%%%%%%%
      toggle_analysis_value("apply_bps");

      % Update all figures (they will all be affected)
      for n = 1:length(plot_names)
        figures.(plot_names{n}).drawfcn();
      end % for

    case 'i'
      %%%%%%%%%%%%%%%%%%%
      % Invert colormap %
      %%%%%%%%%%%%%%%%%%%
      plots.(plot_name).cinv = ~plots.(plot_name).cinv;
      set_unsaved_changes(true);
      figures.(plot_name).drawfcn();

    case 'k'
      %%%%%%%%%%%%%%%%%%%%%%%
      % Select profile mask %
      %%%%%%%%%%%%%%%%%%%%%%%

      % Only do anything if the profile plot is open
      if (isempty(figures.profile.ax_handle))
        return;
      end % if

      analysed.profile_title = "Choose low phase for start of mask";
      set_title("profile");
      analysed.profilemask_click = 1;
      set(figures.profile.fig_handle, "windowbuttondownfcn", @select_profilemask_click);

    case 'K'
      %%%%%%%%%%%%%%%%%%%%%%
      % Clear profile mask %
      %%%%%%%%%%%%%%%%%%%%%%

      % Only do anything if the profile plot is open
      if (isempty(figures.profile.ax_handle))
        return;
      end % if

      if (~isempty(analysis.profile_mask))
        analysis.profile_mask = [];
        set_unsaved_changes(true);
        % Redraw all the figures
        for n = 1:length(plot_names)
          figures.(plot_names{n}).drawfcn();
        end % for
      end % if

    case 'l'
      %%%%%%%%%%%%%%%%%%%%%%%
      % Toggle log plotting %
      %%%%%%%%%%%%%%%%%%%%%%%
      if (~figures.(plot_name).isfft)
        return % (Nothing to be done for non-FFT plots)
      end % if

      % Toggle the islog value
      plots.(plot_name).islog = ~plots.(plot_name).islog;
      set_unsaved_changes(true);

      % Get current axis limits depending on whether it's a 1D or 2D plot
      switch figures.(plot_name).dims
        case 1 % 1D
          ax = plots.(plot_name).axis(3:4);
        case 2 % 2D
          ax = plots.(plot_name).cax;
      end % switch

      % If values were negative, enforce positive
      if (plots.(plot_name).islog)
        if (ax(1) <= 0)
          ax(1) = max([min(abs(analysed.spectrum_vals)), eps]);
        end % if
        if (ax(2) <= 0)
          ax(2) = max([max(abs(analysed.spectrum_vals)), eps]);
        end % if
      end % if

      % Have to change c-axis manually
      if (figures.(plot_name).dims == 2)
        if (plots.(plot_name).islog)
          plots.(plot_name).cax = log10(ax);
        else
          plots.(plot_name).cax = 10.^ax;
        end % if
      else
        plots.(plot_name).axis(3:4) = ax;
      end % if

      % Replot just this plot
      figures.(plot_name).drawfcn();

    case 'm'
      %%%%%%%%%%%%%%%%%%%%%%%%%
      % Toggle Hamming window %
      %%%%%%%%%%%%%%%%%%%%%%%%%
      toggle_analysis_value("apply_hamming");
      figures.fft.drawfcn();

    case 'n'
      %%%%%%%%%%%%%%%%%%%%%%%%%
      % Toggle Hanning window %
      %%%%%%%%%%%%%%%%%%%%%%%%%
      toggle_analysis_value("apply_hanning");
      figures.fft.drawfcn();

    case 'o'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Change period by FFT click %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (~isempty(figures.fft.fig_handle))
        analysed.fft_title = get(get(figures.fft.ax_handle, "title"), "string"); % Save the old title
        title(figures.fft.ax_handle, "Click on a harmonic of the desired frequency\n(right button to cancel)");
        set(figures.fft.fig_handle, "windowbuttondownfcn", @select_pulsar_click);
      end % if

    case 'p'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Change folding period manually %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      cstr = inputdlg({"Enter period in seconds:"}, "Period", 1, {num2str(analysis.period, 15)});
      if (~isempty(cstr)) % = OK was pushed
        set_period(str2num(cstr{1}));
      end % if

    case 'P'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Read in period (1/F0) from ephemeris  (BWM) %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      % Open up an Open File dialog box
      [parfile, parpath] = uigetfile({"*.par", "Ephemeris file"});

      if (~strcmp(parpath, "0")) % then they actually selected a file

        parvals = read_par([parpath, parfile]);

        if isfield(parvals,'P0')
          set_period(str2num(parvals.P0))
        elseif isfield(parvals,'F0')
          set_period(1/str2num(parvals.F0))
        else
          warndlg("No rotation frequency (F0) or period (P0) found in ephemeris.")
        end % if

      end % if

      %if isfield(parvals,'F1')
      %  set_pdot(1/str2num(parvals.F1'))
      %else
      %  set_pdot(0)   
      %end % if

    case 'q'
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Toggle quantisation when adding filters %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (isfield(analysed, "filter") && strcmp(plot_name, "tdfs")) % This only happens when in "create filter" mode
        switch analysed.filter.quantised
          case "on"
            analysed.filter.quantised = "off";
          case "off"
            analysed.filter.quantised = "on";
        end % switch
        analysed.tdfs_title = ["Click on the centre of the new (", analysed.filter.direction, ...
                               ") filter\nPress 'd' to toggle direction filter\n", ...
                               "Press 'q' to toggle quantisation (currently ", ...
                               analysed.filter.quantised, ")\n", ...
                               "Right click to cancel"];
        set_title(plot_name);
      end % if

    case 'v'
      %%%%%%%%%%%%%%%%%%%%%%%
      % Toggle only_visible %
      %%%%%%%%%%%%%%%%%%%%%%%
      toggle_analysis_value("only_visible");
      figures.fft.drawfcn();
      % Does profile plot need to be redrawn?
      figures.pulsestack.drawfcn();
      figures.tdfs.drawfcn();
      figures.modenv.drawfcn();

    case 'V'
      %%%%%%%%%%%%%%%%%%%%%%%
      % Toggle only_visible %
      %%%%%%%%%%%%%%%%%%%%%%%
      toggle_analysis_value("only_visible_stack");
      figures.tdfs.drawfcn();
      figures.modenv.drawfcn();

    case 's'
      %%%%%%%%%%%%%%%%%%%%%%%%
      % Change sampling rate %
      %%%%%%%%%%%%%%%%%%%%%%%%

      % Read in the new value via a dialog box
      cstr = inputdlg({["Please enter the sampling rate (", data.frequnits, "):"]}, ...
                       "Set sampling rate", 1, {num2str(data.samplingrate,15)});

      if (~isempty(cstr))
        new_samplingrate = str2num(cstr{1});

        % If the new value is different from the old
        if (data.samplingrate ~= new_samplingrate)
          % Record the scale factor by which the sampling rate has changed
          scale = new_samplingrate / data.samplingrate;

          % Update values
          if (~isempty(analysis.breakpoints))
            analysis.breakpoints *= data.samplingrate / new_samplingrate;
          end % if
          data.samplingrate = new_samplingrate;
          analysed.dt = 1/data.samplingrate;                 % The time between adjacent samples
          analysed.t  = [0:(analysed.N-1)]' * analysed.dt;   % The time axis

          % Rescale x axis of affected plots
          if (~isempty(plots.timeseries.axis))
            plots.timeseries.axis(1:2) /= scale;
          end % if
          if (~isempty(plots.fft.axis))
            plots.fft.axis(1:2)        *= scale;
          end % if

          % Update all figures (they will all be affected)
          for n = 1:length(plot_names)
            figures.(plot_names{n}).drawfcn();
          end % for

          % Set unsaved changes flag
          set_unsaved_changes(true);

        end % if
      end % if

    case 'z'
      %%%%%%%%%%%%%%%%%%%
      % Toggle zeromean %
      %%%%%%%%%%%%%%%%%%%
      toggle_analysis_value("zeromean");
      figures.fft.drawfcn();

  end % switch

end % function

function this_plot_name = fighandle2plotname(f)

  global figures;

  % Get the plot_name of the figure window in which the key was pressed
  this_plot_name = [];
  plot_names = fieldnames(figures);
  for n = 1:length(plot_names)
    plot_name = plot_names{n};
    if (f == figures.(plot_name).fig_handle)
      this_plot_name = plot_name;
      break;
    end % if
  end % for

end % function

function set_analysis_value(name, value)

  global analysis;

  if (~isfield(analysis, name))
    warndlg(["Creating new field \"", name, "\""]);
    analysis.(name) = [];
  end % if

  if (~isequal(analysis.(name), value))
    analysis.(name) = value;
    set_unsaved_changes(true);
  end % if

end % function

function toggle_analysis_value(name)

  global analysis;

  set_analysis_value(name, ~analysis.(name));

end % function

function collect_breakpoint_clicks(src, button)

  global figures;
  global analysis;
  global analysed;

  action = [];

  switch button
    case 1 % left mouse button
      action = @add_breakpoint;
    case 3 % right mouse button
      if (length(analysis.breakpoints) > 0)
        action = @remove_breakpoint;
      end % if
  end % switch

  if (~isempty(action))
    point = get(figures.timeseries.ax_handle, "currentpoint");
    action(point(1));
    flatten();
    figures.timeseries.drawfcn();
  end % if

end % function

function add_breakpoint(x)
  global analysis;

  analysis.breakpoints = union(analysis.breakpoints, x); % <-- add a breakpoint
end % function

function remove_breakpoint(x)
  global analysis;

  [dummy, nearest_idx] = min(abs(analysis.breakpoints - x));
  analysis.breakpoints = setdiff(analysis.breakpoints, analysis.breakpoints(nearest_idx)); % <-- remove a breakpoint
end % function

function select_pulsar_click(src, button)

  global figures;
  global analysed;

  if (button == 1) % Left button = convert clicked frequency to period
    cstr = inputdlg({"Harmonic number of selected point:"}, "Harmonic", 1, {"1"});
    if (~isempty(cstr))

      % Get where the user clicked
      point = get(figures.fft.ax_handle, "currentpoint");
      nharm = str2num(cstr{1});
      newperiod = nharm / point(1,1);
      set_period(newperiod);

    end % if
  end % if

  if (any(button == [1,3])) % Reset callback function as long as either the left or right buttons were pushed
    title(figures.fft.ax_handle, analysed.fft_title);
    set(figures.fft.fig_handle, "windowbuttondownfcn", []);
    analysed = rmfield(analysed, "fft_title");
  end % if

end % function

function set_period(newperiod)

  global figures;

  set_analysis_value("period", newperiod);

  % Update all figures (they will all be affected)
  plot_names = fieldnames(figures);
  for n = 1:length(plot_names)
    figures.(plot_names{n}).drawfcn();
  end % for

  % Update the title on the profile plot
  set_title("profile");

end % function

function select_profilemask_click(src, button)

  global figures;

  global analysis;
  global analysed;

  f = figures.profile.fig_handle;
  a = figures.profile.ax_handle;

  switch button
    case 1 % Left mouse button
      point = get(a, "currentpoint"); % Get click coordinates
      switch analysed.profilemask_click
        case 1 % on first click
          analysed.profile_mask = point(1,1); % Save clicked point
          analysed.profile_title = "Choose high phase for end of mask";
          set_title("profile");
          analysed.profilemask_click = 2;
        case 2 % on second click
          set_analysis_value("profile_mask", [analysed.profile_mask, point(1,1)]);
          analysed = rmfield(analysed, {"profile_title", "profile_mask", "profilemask_click"});
          % Redraw all the figures
          plot_names = fieldnames(figures);
          for n = 1:length(plot_names)
            figures.(plot_names{n}).drawfcn();
          end % for
          set(f, "windowbuttondownfcn", []);
      end % switch

    case 3 % Right mouse button

  end % switch

end % function

%%%%%%%%%%%%%%%%%%%%%%%%
% PAN / ZOOM FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%

function panzoom_down(src, button, plot_name)

  global analysed;
  global figures;

  % Only do anything if the left or right mouse buttons were clicked
  if (~any(button == [1,3]))
    return
  end % if

  % Only do anything if a pan-zoom-able axes was clicked
  a = gca();
  cont = false; % Flag for whether to continue or not
  plot_names = fieldnames(figures);
  for n = 1:length(plot_names)
    if (any(a == figures.(plot_names{n}).ax_handle))
      cont = true;
    end % if
  end % for

  if (~cont)
    return
  end % if

  % Get the point that was clicked in axis coordinates, and which button was used
  analysed.panzoom.point1 = fig2ax_coords(plot_name);
  analysed.panzoom.button = button;

  % Set the motion and up callback functions
  set(src, "windowbuttonmotionfcn", {@panzoom_motion, plot_name});
  set(src, "windowbuttonupfcn",     {@panzoom_up,     plot_name});

end % function

function panzoom_motion(src, button, plot_name)

  global figures;
  global plots;
  global analysed;

  % Get the current cursor position
  analysed.panzoom.point2 = fig2ax_coords(plot_name);

  % Extra processing (namely, moving the axes), if we're in "pan" mode
  if (analysed.panzoom.button == 1) % = left button = pan
    % If y-axis is logscale, calculate things accordingly
    islog = false;
    if (isfield(plots.(plot_name), "islog"))
      if (plots.(plot_name).islog && (figures.(plot_name).dims == 1))
        % Assumption: only y-axes that fit these conditions are allowed to be logscale
        islog = true;
      end % if
    end % if

    % Get subplot number
    s = find(gca == figures.(plot_name).ax_handle);

    % (Linear) change in position
    delta_pos = analysed.panzoom.point2 - analysed.panzoom.point1;
    plots.(plot_name).axis(s,1:2) -= delta_pos(1); % Change the x axis in any case

    if (~islog)
      plots.(plot_name).axis(s,3:4) -= delta_pos(2);
    else
      delta_pos_log = analysed.panzoom.point2 ./ analysed.panzoom.point1;
      plots.(plot_name).axis(s,3:4) ./= delta_pos_log(2);
    end % if

    if (strcmp(plot_name, "modenv"))
      % Get "corresponding" plot, i.e. with matching x-axis
      s2 = mod(s+1,4)+1; % 1->3, 3->1, 2->4, 4->2
      % Make x axes equal
      plots.(plot_name).axis(s2,1:2) = plots.(plot_name).axis(s,1:2);
    end % if

    set_unsaved_changes(true);

    % Don't replot the figure, just set the axes
    for n = 1:figures.(plot_name).nplots
      axis(figures.(plot_name).ax_handle(n), plots.(plot_name).axis(n,:));
    end % for

  end % if

end % function

function panzoom_up(src, button, plot_name)

  global figures;
  global plots;
  global analysis;
  global analysed;

  % Panning takes care of itself, but when zooming...
  if (analysed.panzoom.button == 3) % = right button = zoom

    if (~isfield(analysed.panzoom, "point2")) % then no mouse drag occurred
      % Clean up and get outta here
      % Clear relevant "analysed" variables
      analysed = rmfield(analysed, "panzoom");

      % Set the motion and up callback functions
      set(src, "windowbuttonmotionfcn", []);
      set(src, "windowbuttonupfcn",     []);

      return
    end % if

    % Get subplot number
    s = find(gca == figures.(plot_name).ax_handle);

    xmin = min([analysed.panzoom.point1(1), analysed.panzoom.point2(1)]);
    xmax = max([analysed.panzoom.point1(1), analysed.panzoom.point2(1)]);
    ymin = min([analysed.panzoom.point1(2), analysed.panzoom.point2(2)]);
    ymax = max([analysed.panzoom.point1(2), analysed.panzoom.point2(2)]);

    % Only zoom if zoom window has non-zero height and width
    if ((xmax > xmin) && (ymax > ymin))
      plots.(plot_name).axis(s,:) = [xmin, xmax, ymin, ymax];
      set_unsaved_changes(true);
    end % if

    if (strcmp(plot_name, "modenv"))
      % Get "corresponding" plot, i.e. with matching x-axis
      s2 = mod(s+1,4)+1; % 1->3, 3->1, 2->4, 4->2
      % Make x axes equal
      plots.(plot_name).axis(s2,1:2) = plots.(plot_name).axis(s,1:2);
    end % if

    % Don't replot the figure, just set the axes
    for n = 1:figures.(plot_name).nplots
      axis(figures.(plot_name).ax_handle(n), plots.(plot_name).axis(n,:));
    end % for

  end % if

  % Clear relevant "analysed" variables
  analysed = rmfield(analysed, "panzoom");

  % Set the motion and up callback functions
  set(src, "windowbuttonmotionfcn", []);
  set(src, "windowbuttonupfcn",     []);

  % If they've moved the timeseries plot and only_visible is turned on,
  % then replot everything.
  if (strcmp(plot_name, "timeseries") && analysis.only_visible)
    plot_names = fieldnames(figures);
    for n = 1:length(plot_names)
      figures.(plot_names{n}).drawfcn();
    end % for
  end % if

  % If they've moved the pulsestack plot and only_visible_stack is turned on,
  % then replot 2DFS and modulation envelope plots
  if (strcmp(plot_name, "pulsestack") && analysis.only_visible_stack)
    figures.tdfs.drawfcn();
    figures.modenv.drawfcn();
  end % if

end % function

function point = fig2ax_coords(plot_name)

  global figures;
  global plots;

  f = figures.(plot_name).fig_handle;
  a = gca; % <-- I have to use this instead of my own figures.(plot_name).ax_handle in order
           %     to make it work for figures with more than one plot, namely, "modenv".

  % Get subplot number (will equal 1 for figures with only one plot)
  s = find(a == figures.(plot_name).ax_handle);

  % If either figure or axes handles don't exist, don't bother
  if (isempty(f) || isempty(a))
    return
  end % if

  % If y-axis is logscale, calculate things accordingly
  islog = false;
  if (isfield(plots.(plot_name), "islog"))
    if (plots.(plot_name).islog && (figures.(plot_name).dims == 1))
      % Assumption: only y-axes that fit these conditions are allowed to be logscale
      islog = true;
    end % if
  end % if

  point   = get(f, "currentpoint");
  f_pos   = get(f, "position");
  f_width = f_pos(3:4);
  a_pos   = get(a, "position"); % values are in figure window units
  a_axis  = plots.(plot_name).axis(s,:);
  if (islog)
    a_axis(3:4) = log(a_axis(3:4));
  end % if
  a_width = diff(reshape(a_axis,2,2));
  a_orig  = a_axis([1,3]);
  point   = (point./f_width - a_pos(1:2)) ./ a_pos(3:4).*a_width + a_orig;
  if (islog)
    point(2) = exp(point(2));
  end % if

end % function

function add_filter(src, button)

  global figures;
  global analysis;
  global analysed;

  switch button
    case 1 % Left button clicked
      % Get the last clicked position
      a = figures.tdfs.ax_handle;
      pos = get(a, "currentpoint");
      pos = pos(1,1:2);

      % Quantise, if necessary
      if (strcmp(analysed.filter.quantised, "on"))
        pos = quantise_tdfs_pos(pos);
      end % if

      switch analysed.filter.click_no
        case 1 % On first click, choose centre of filter
          analysed.filter.pos = pos;
          analysed.tdfs_title = ["Click on the edge of the new (", analysed.filter.direction, ...
                                 ") filter\nPress 'd' to toggle direction filter\n", ...
                                 "Press 'q' to toggle quantisation (currently ", ...
                                 analysed.filter.quantised, ")\n", ...
                                 "Right click to cancel"];
          set_title("tdfs");
          analysed.filter.click_no = 2;
        case 2 % On second click, choose edge of filter
          % Calculate filter parameters
          switch analysed.filter.direction
            case "horizontal"
              centre = analysed.filter.pos(2);
              width  = abs(pos(2) - analysed.filter.pos(2));
              dir    = 0;
            case "vertical"
              centre = analysed.filter.pos(1);
              width  = abs(pos(1) - analysed.filter.pos(1));
              dir    = 1;
          end % switch

          % Add filter to list
          if (width > 0)
            analysis.filters = [analysis.filters; [centre, width, dir]];
            set_unsaved_changes(true);
          else
            errordlg("Cannot choose filter with zero width");
          end % if

          % Reset title
          analysed = rmfield(analysed, "tdfs_title");

          % Reset mouse click callback
          set(figures.tdfs.fig_handle, "windowbuttondownfcn", {@panzoom_down, "tdfs"});

          % Clean up "analysed" structure
          analysed = rmfield(analysed, "filter");

          % Redraw plots
          figures.tdfs.drawfcn();
          figures.modenv.drawfcn();
      end % switch

    case 3 % Right button clicked

      % Clean up and get outta here
      analysed = rmfield(analysed, {"filter", "tdfs_title"});
      figures.tdfs.drawfcn();
      set(figures.tdfs.fig_handle, "windowbuttondownfcn", {@panzoom_down, "tdfs"});

  end % switch

end % function

function remove_filter(src, button)

  global figures;
  global analysis;
  global analysed;

  switch button
    case 1 % Left button
      % Get the last clicked position
      a = figures.tdfs.ax_handle;
      pos = get(a, "currentpoint");
      pos = pos(1,1:2);

      % Loop through filters and delete any containing the clicked region
      to_be_kept = [];
      for n = 1:rows(analysis.filters)
        centre    = analysis.filters(n,1);
        width     = analysis.filters(n,2);
        direction = analysis.filters(n,3);

        switch direction
          case 0 % horizontal
            pos_coord = pos(2);
          case 1 % vertical
            pos_coord = pos(1);
        end % switch

        % Check if user clicked OUTside filter
        if (abs(pos_coord - centre) > width)
          to_be_kept = [to_be_kept, n];
        end % if
      end % for

      % Only do anything if they successfully clicked inside at least one filter
      if (length(to_be_kept) == rows(analysis.filters))
        return
      end % if

      analysis.filters = analysis.filters(to_be_kept,:);
      set_unsaved_changes(true);

      % Reset title
      analysed = rmfield(analysed, "tdfs_title");

      % Clean up "analysed" structure
      analysed = rmfield(analysed, "filter");

      % Redraw plots
      figures.tdfs.drawfcn();
      figures.modenv.drawfcn();
 
      % Reset mouse click callback
      set(figures.tdfs.fig_handle, "windowbuttondownfcn", {@panzoom_down, "tdfs"});

    case 3 % Right button

      % Clean up "analysed" structure
      analysed = rmfield(analysed, "filter");

      % Reset title
      analysed = rmfield(analysed, "tdfs_title");
      figures.tdfs.drawfcn();

      % Reset mouse click callback
      set(figures.tdfs.fig_handle, "windowbuttondownfcn", {@panzoom_down, "tdfs"});

  end % switch

end % function

function shift_DC(src, button)

  global figures;
  global analysis;
  global analysed;

  switch button
    case 1 % Left button

      % Get the last clicked position
      a = figures.tdfs.ax_handle;
      pos = get(a, "currentpoint");
      analysis.shift_DC = quantise_tdfs_pos(pos(1,1:2));

      % Reset title
      analysed = rmfield(analysed, "tdfs_title");

      % Reset mouse click callback
      set(figures.tdfs.fig_handle, "windowbuttondownfcn", {@panzoom_down, "tdfs"});

      % Redraw figure
      figures.tdfs.drawfcn();
      figures.modenv.drawfcn();

    case 3 % Right button

      % Reset title
      analysed = rmfield(analysed, "tdfs_title");
      set_title("tdfs");

      % Reset mouse click callback
      set(figures.tdfs.fig_handle, "windowbuttondownfcn", {@panzoom_down, "tdfs"});

  end % switch

end % function

function new_pos = quantise_tdfs_pos(pos)

  global analysed;

  dxy = [analysed.tdfs.dx, analysed.tdfs.dy];
  new_pos = round(pos./dxy) .* dxy;

end % function
