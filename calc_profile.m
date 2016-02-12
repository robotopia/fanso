function calc_profile()

  global data;
  global analysis;
  global analysed;

  if (analysis.apply_bps)
    to_be_folded = analysed.flattened;
  else
    to_be_folded = data.timeseries;
  end % if

  % Calculate the values for the profile abscissa (=phase)
  N = length(to_be_folded);
  t = analysed.t + 0.5*analysed.dt; % +0.5 is to avoid some of the more common computer precision errors
  phase = mod(t, analysis.period) / analysis.period;
  accum_subs = floor(phase * analysis.nprofile_bins) + 1;
  analysed.profile = accumarray(accum_subs, to_be_folded, [analysis.nprofile_bins,1], @mean);

end % function
