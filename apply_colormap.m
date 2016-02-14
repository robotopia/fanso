function apply_colormap(plot_name)

  global figures;
  global plots;

  if (figures.(plot_name).dims ~= 2) % then not a 2D plot
    return
  end % if

  % Get axis handle
  a = figures.(plot_name).ax_handle;

  colormap_name = colormap("list"){plots.(plot_name).cmap};
  to_be_inverted = plots.(plot_name).cinv;
  cmap = colormap(a, colormap_name);
  if (to_be_inverted)
    colormap(a, flipud(cmap));
  end % if

end % function
