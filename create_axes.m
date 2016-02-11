function a = create_axes(plot_name)

  global figures;

  % Only create axes if the figure window already exists
  f = figures.(plot_name).fig_handle;
  if (isempty(f))
    a = [];
    return
  end % if
  figure(f);

  % Set/Create axes
  a = figures.(plot_name).ax_handle;
  if (isempty(a))
    a = axes();
    figures.(plot_name).ax_handle = a;

  end % if

end % function
