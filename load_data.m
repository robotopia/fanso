function load_data(loadpath, loadfile)

  try
    % Load file contents
    load("-binary", [loadpath, loadfile]);

    % Store file name + path in global variables
    global filepath;
    global filename;

    filename = loadfile;
    filepath = loadpath;

  catch
    errordlg({["Unable to load file \"", loadpath, "\", loadfile, "\""], lasterr()});
  end % try_catch

end % function
