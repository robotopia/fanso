function stackpulses(applytransforms)

  global data;
  global analysis;
  global analysed;

  if (~applytransforms)
    % Still zero-pad (so that we can form a rectangle matrix)
    np = analysis.period / analysed.dt; % Number of bins per pulse
    n_extra_zeros = round(ceil(analysed.N/np)*np) - analysed.N;
    to_be_stacked = [analysed.flattened; zeros(n_extra_zeros,1)];
  else
    apply_transforms();
    to_be_stacked = analysed.transformed;
  end % if

  % Prepare the grid of values
  N         = length(to_be_stacked);
  t         = ([1:N]' - 0.5) * analysed.dt; % -0.5 is to avoid some of the more common computer precision errors
  pulse_no  = floor(t/analysis.period) + 1;
  phase_bin = floor(mod(t,analysis.period)/analysed.dt) + 1; % The "floor" here is problematic. It means some of the pulses
                                                             % may be shifted by up to dt/2. I see no other way around this.
  accum_subs = [pulse_no, phase_bin];
  stacked    = accumarray(accum_subs, to_be_stacked);

  if (~applytransforms)
    analysed.stacked.all = stacked;
  else
    analysed.stacked.transformed = stacked;
  end % if

end % function
