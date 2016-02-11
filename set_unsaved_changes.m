function set_unsaved_changes(newval)

  global unsaved_changes
  global figures
  global filename

  unsaved_changes = newval;

  h = figures.timeseries.fig_handle;
  figure(h);

  fig_name = ["timeseries"];

  if (filename)
    fig_name = [fig_name, ": ", filename];
  end % if

  if (unsaved_changes)
    % Add a "*" if there are unsaved changes
    fig_name = [fig_name, "*"];
  end

  figure(h, "name", fig_name);

end % function

