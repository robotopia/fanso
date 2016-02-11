function load_fan(loadfile)

  global filepath;
  global filename;

  global figures;

  % Ask about unsaved changes
  if (~offer_to_save())
    return
  end % if

  % Make sure loadfile is a string before testing its contents
  if (~ischar(loadfile))
    loadfile = "";
  end % if

  if (exist(loadfile, "file") == 2) % We already have an existing file
    loadpath = "";
  else % Open up an Open File dialog box
    [loadfile, loadpath] = uigetfile({"*.fan", "FANSO file"});
  end

  if (~strcmp(loadpath, "0")) % then they actually selected a file
    try
      % Load file contents
      load("-binary", [loadpath, loadfile]);
    catch
      errordlg({["Unable to load file \"", loadfile, "\""], lasterr()});
    end % try_catch

    % Store file name + path in global variables
    filename = loadfile;
    filepath = loadpath;

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
