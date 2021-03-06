function plot_modenv()

  global figures;

  global data;
  global plots;
  global analysis;
  global analysed;

  plot_name = "modenv";

  a = create_axes(plot_name);
  if (isempty(a))
    return
  end % if

  % If 2DFS hasn't yet been calculated, plot tdfs first
  if (isempty(figures.tdfs.fig_handle))
    create_figure("tdfs");
    figures.tdfs.drawfcn();
    figure(figures.(plot_name).fig_handle);
  end % if

  % This is a special plot, which needs special handling
  a1 = subplot(2,2,1);
  a2 = subplot(2,2,2);
  a3 = subplot(2,2,3);
  a4 = subplot(2,2,4);
  figures.(plot_name).ax_handle = [a1, a2, a3, a4];
  
  ax = plots.(plot_name).axis; % <-- This will have 4 rows in it

  % Unshift filtered plot, where the shift is due to both centering and user-defined DC-shift
  xshift = -analysed.tdfs.xshift;
  yshift = -analysed.tdfs.yshift;

  if (~isempty(analysis.shift_DC))
    % Convert shift_DC amounts to pixel units
    dy = 1/rows(analysed.tdfs.filtered);
    shift_DC = round(analysis.shift_DC ./ [1, dy]);
    xshift -= shift_DC(1);
    yshift -= shift_DC(2);
  end % if

  analysed.tdfs.unshifted = shift(analysed.tdfs.filtered,  xshift, 2);
  analysed.tdfs.unshifted = shift(analysed.tdfs.unshifted, yshift, 1);

  % Now calculate the inverse 2DFS (=m_{phi,t} in E&S(2002))
  M = ifft2(2*analysed.tdfs.unshifted); % For explanation of extra factor of 2, see E&S (2002)

  % Do singular-value decomposition on m_{phi,t} = m_phi(phi) * m_t(t)
  [U,S,V] = svd(M);
  if (~isfield(analysed, "Sn"))
    analysed.Sn = 1;
    analysed.SN = min(size(S));
  end % if
  Sn = analysed.Sn;
  m_t     = U(:,Sn) * sqrt(S(Sn,Sn));
  m_phi   = V(:,Sn) * sqrt(S(Sn,Sn));

  % Put U,S,V somewhere globally accessible
  analysed.modenv.U = U;
  analysed.modenv.S = S;
  analysed.modenv.V = V;

  % Calculate various x-axes
  if (analysis.only_visible_stack)
    pax = plots.pulsestack.axis;
    tmin = max([floor(pax(3)), 1]);
    pmin_idx = max([floor(pax(1) * columns(analysed.stacked.all)), 0]) + 1;
fflush(stdout);
  elseif (analysis.only_visible)
    tax = plots.timeseries.axis;
    tmin = max([floor(tax(1)/analysis.period),0]) + 1;
    pmin_idx = 1;
  else
    tmin = 1;
    pmin_idx = 1;
  end % if

  tmax     = tmin     + length(m_t)   - 1;
  pmax_idx = pmin_idx + length(m_phi) - 1;

  pulse_no  = [tmin:tmax];
  longitude = [pmin_idx:pmax_idx] / columns(analysed.stacked.all) * 360;

  y1 = abs(m_t);
  y2 = abs(m_phi);
  y3 = arg(m_t) * 180/pi;
  y4 = arg(m_phi) * 180/pi;

  plots.(plot_name).autoscale = [min(pulse_no),  max(pulse_no),  min(y1), max(y1);
                                 min(longitude), max(longitude), min(y2), max(y2);
                                 min(pulse_no),  max(pulse_no),  min(y3), max(y3);
                                 min(longitude), max(longitude), min(y4), max(y4)];

  % Pulse no. vs Amplitude
  plot(a1, pulse_no, y1);
  xlabel(a1, "Pulse Number");
  ylabel(a1, "Amplitude");

  % Longitude vs Amplitude
  plot(a2, longitude, y2);
  xlabel(a2, "Pulse longitude");
  ylabel(a2, "Amplitude");

  % Pulse no. vs Phase
  plot(a3, pulse_no, y3, 'x');
  xlabel(a3, "Pulse Number");
  ylabel(a3, "Phase");

  % Longitude vs Phase
  plot(a4, longitude, y4, 'x');
  xlabel(a4, "Pulse longitude");
  ylabel(a4, "Phase");

  % Calculate the "signal-to-noise" of the result
  analysed.SNR = S(Sn,Sn) / (trace(S) - S(Sn,Sn)); % <-- Not sure if this is the correct way to do it...
  set_title(plot_name);

  % Get/Set axis limits
  if (~isempty(ax))
    axis(a1, ax(1,:));
    axis(a2, ax(2,:));
    axis(a3, ax(3,:));
    axis(a4, ax(4,:));
  else
    plots.(plot_name).axis = [axis(a1); axis(a2); axis(a3); axis(a4)];
  end % if

end % function
