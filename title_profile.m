function title_profile()

  global figures;

  global data;
  global analysis;
  global analysed;

  a = figures.profile.ax_handle;

  if (~isempty(a))
    if (isfield(analysed, "profile_title"))
      title(a, analysed.profile_title);
    else
      title(a, sprintf("Period = %.12f %s;    No. of bins = %d", analysis.period, data.timeunits, analysis.nprofile_bins));
    end % if
  end % if

end % function
