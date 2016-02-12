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

    % Set callback function for when a key is pressed
    if (~isempty(data))
      set(figures.(plot_name).fig_handle, "keypressfcn", @keypressfcn);
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

  % Set enables on menu items
  global m_data_save m_data_saveas m_data_export;
  set([m_data_save, m_data_saveas, m_data_export], "enable", "off");

  % Set callback function for when a key is pressed
  set(figures.timeseries.fig_handle, "keypressfcn", []);
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
function saveas_fan()

  global filepath;
  global filename;

  % Open up a Save File dialog box
  if (filepath && filename)
    [savefile, savepath] = uiputfile([filepath, filename]);
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
  end % if

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function import_timeseries()

  global figures;

  global data;
  global plots;
  global analysis;

  global filename;
  global filepath;

  % Ask about unsaved changes
  if (~offer_to_save())
    return
  end % if

  % Get file from Open File dialog box
  [loadfile, loadpath] = uigetfile();

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

    % Reset all the other variables
    % "plots" structure
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
    analysis.profile_mask   = []; % A pair of phases that define a region of phases to be ignored in the breakpoint linear fits
    analysis.period         = []; % The folding period
    analysis.P2hat          = []; % The measured longitudinal "time" between subpulses, P2
    analysis.P3hat          = []; % The measured "time" between subpulses at the same phase, P3
    analysis.nprofile_bins  = []; % The number of bins to be used for folding output
    analysis.zeromean       = 0;  % 0 = Do nothing;                1 = Zero mean before applying FFT
    analysis.zeropad        = 0;  % 0 = Do nothing;                1 = Zero-pad to a nearly-whole number of periods
    analysis.only_visible   = 0;  % 0 = Analyse entire timeseries; 1 = Analyse only visible timeseries
    analysis.apply_hamming  = 0;  % 0 = Do nothing;                1 = Apply Hamming window
    analysis.apply_hanning  = 0;  % 0 = Do nothing;                1 = Apply Hanning window
    analysis.apply_bps      = 0;  % 0 = Do nothing;                1 = Apply breakpoints (i.e. "flatten" timeseries)
    analysis.breakpoints    = []; % The x-coordinates of the breakpoints
    analysis.filters        = []; % Horizontal and vertical filters (E&S 2002 - style) used on the 2DFS
    analysis.shift_DC       = []; % Horizontal and vertical displacement of the origin in the 2DFS

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

    % Clear the (saved) filename and path variables
    filename = [];
    filepath = [];

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

  global filepath;
  global filename;
% SOMETHING WRONG IN THIS FunCTION.
% load file
% make change
% choose New
% it goes to "save as" instead of "save"

  global unsaved_changes;
  cont = true;

  if (unsaved_changes)
    dlg_answer = questdlg("Would you like to save your changes?", "Save changes?");

    switch dlg_answer
      case "Yes"
        if (filename && filepath)
          save_fan();
        else
          saveas_fan();
        end
      case "Cancel"
        cont = false;
    end % switch
  end % if

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function close_figure(src, data, plot_name)

  global figures

  % Save current views
  save_views();

  %if (any(strcmp(evt.Key, {"f1", "f2", "f3", "f4", "f5", "f6", "f7"})))
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

  global figures

  global data
  global plots
  global analysis
  global analysed

  % Save current views
  save_views();

  % Get the plot_name of the figure window in which the key was pressed
  this_plot_name = [];
  plot_names = fieldnames(figures);
  for n = 1:length(plot_names)
    plot_name = plot_names{n};
    if (src == figures.(plot_name).fig_handle)
      this_plot_name = plot_name;
      break;
    end % if
  end % for

  if (isempty(this_plot_name))
    error("Unknown callback figure handle in keypressfcn");
  end % if

  %if (any(strcmp(evt.Key, {"f1", "f2", "f3", "f4", "f5", "f6", "f7"})))
  if (any(strcmp(evt.Key, {"f1", "f2", "f3"})))
    plot_no    = str2num(evt.Key(2));
    plot_name  = plot_names{plot_no};
    if (isempty(figures.(plot_name).fig_handle))
      create_figure(plot_name);
      figures.(plot_name).drawfcn();
    end % if
  end % switch

  switch evt.Character
    case 'h'
      %%%%%%%%%%%%%%%%
      % Display help %
      %%%%%%%%%%%%%%%%

      msgbox({"F1 = plot timeseries",
              "F2 = plot FFT",
              "F3 = plot profile",
              "F4 = plot HRFS",
              "F5 = plot waterfall",
              "F6 = plot 2DFS",
              "F7 = plot modulation envelopes",
              "0 = toggle zero-padding",
              "^ = toggle show peaks",
              "b = change number of profile bins",
              "B = toggle breakpoint edit mode",
              "f = flatten timeseries",
              "h = display this help",
              "l = toggle logarithmic plot for current figure",
              "m = toggle Hamming window",
              "n = toggle Hanning window",
              "o = change folding period by clicking on FFT",
              "p = change folding period manually",
              "P = change folding period by selecting pulsar",
              "s = change sampling rate",
              "v = toggle only analyse visible part of timeseries",
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

    case '^'
      %%%%%%%%%%%%%%%%%%%%%
      % Toggle show peaks %
      %%%%%%%%%%%%%%%%%%%%%
      toggle_analysis_value("show_peaks");
      figures.tdfs.drawfcn();

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
        title(figures.timeseries.ax_handle, "");

        % Reset callback functions
        set(figures.timeseries.fig_handle, "windowbuttondownfcn", []);

      else % i.e. it is NOT in edit mode when b was pressed

        analysed.old_apply_bps   = analysis.apply_bps;
        analysed.old_breakpoints = analysis.breakpoints;
        analysed.title           = "Left mouse button = add breakpoint; right = remove breakpoint";

        if (analysis.apply_bps)
          analysis.apply_bps = 0; % i.e. turn it OFF (temporarily)
          figures.timeseries.drawfcn();
        end % if

        title(figures.timeseries.ax_handle, analysed.title);

        set(figures.timeseries.fig_handle, "windowbuttondownfcn", @collect_breakpoint_clicks);

      end % if

      % Flip the mode
      analysed.bp_editmode = ~analysed.bp_editmode;

    case 'f'
      %%%%%%%%%%%%%%%%%%
      % Toggle flatten %
      %%%%%%%%%%%%%%%%%%
      toggle_analysis_value("apply_bps");

      % Update all figures (they will all be affected)
      for n = 1:length(plot_names)
        figures.(plot_names{n}).drawfcn();
      end % for

    case 'l'
      %%%%%%%%%%%%%%%%%%%%%%%
      % Toggle log plotting %
      %%%%%%%%%%%%%%%%%%%%%%%
      if (~figures.(this_plot_name).isfft)
        return % (Nothing to be done for non-FFT plots)
      end % if

      % Toggle the islog value
      plots.(this_plot_name).islog = ~plots.(this_plot_name).islog;
      set_unsaved_changes(true);

      % Get current axis limits depending on whether it's a 1D or 2D plot
      switch figures.(this_plot_name).dims
        case 1 % 1D
          ax = plots.(this_plot_name).axis(3:4);
        case 2 % 2D
          ax = plots.(this_plot_name).cax;
      end % switch

      % If values were negative, enforce positive
      if (plots.(this_plot_name).islog)
        if (ax(1) <= 0)
          ax(1) = max([min(analysed.spectrum_vals), eps]);
        end % if
        if (ax(2) <= 0)
          ax(2) = max([max(analysed.spectrum_vals), eps]);
        end % if
      end % if

      % Have to change c-axis manually
      if (figures.(this_plot_name).dims == 2)
        if (plots.(this_plot_name).islog)
          plots.(this_plot_name).cax = log10(ax);
        else
          plots.(this_plot_name).cax = 10.^ax;
        end % if
      else
        plots.(this_plot_name).axis(3:4) = ax;
      end % if

      % Replot just this plot
      figures.(this_plot_name).drawfcn();

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
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Change folding period by selecting pulsar from list %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      global pulsars

      pulsarnames = fieldnames(pulsars);
      [sel, ok] = listdlg("ListString", pulsarnames, "SelectionMode", "Single", "Name", "Select pulsar");
      if (ok)
        set_period(pulsars.(pulsarnames{sel}).period);
      end % if

    case 'v'
      %%%%%%%%%%%%%%%%%%%%%%%
      % Toggle only_visible %
      %%%%%%%%%%%%%%%%%%%%%%%
      toggle_analysis_value("only_visible");
      figures.fft.drawfcn();
      figures.waterfall.drawfcn();
      figures.tdfs.drawfcn();
      figures.hrfs.drawfcn();
      % Do other plots need to be redrawn?

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

function save_views()

  global figures;
  global plots;

  plot_names = fieldnames(figures);
  for n = 1:length(plot_names)
    plot_name = plot_names{n};
    a = figures.(plot_name).ax_handle;
    if (~isempty(a))
      ax = axis(a);
      if (plots.(plot_name).axis ~= ax)
        plots.(plot_name).axis = ax;
        set_unsaved_changes(true);
      end % if
    end % if
  end % for

end % function

function set_analysis_value(name, value)

  global analysis;

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
    title(figures.timeseries.ax_handle, analysed.title);
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
  title_profile();

end % function
