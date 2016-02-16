function fanso(init_filename)
% -- Function: fanso ()
%     Analyse a timeseries using the Fourier ANalysis Suite for Octave.
%     Requires FLTK.
%
%     Written by Sam McSweeney, 2016, Creative Commons Licence
%     sammy.mcsweeney@gmail.com

  % First, check if there is already an instance of FANSO open
  global __FANSO_INSTANCE__
  if (isempty(__FANSO_INSTANCE__)) % i.e. there is NO instance already open
    __FANSO_INSTANCE__ = 1;
  else
    % Exit gracefully
    disp('FANSO already running');
    return
  end

  % Require that we are using the correct graphics toolkit (FLTK)
  gtk = 'fltk';
  %gtk = 'fltk';
  if (~any(strcmp(available_graphics_toolkits(), gtk)))
    disp(['This function requires ',gtk,' to be installed.'])
    return;
  else
    graphics_toolkit(gtk);
  end % if

  % No file is loaded yet...
  global filename
  filename = [];

  % Setup figure and plot handling system
  global figures;
  figures = struct();

  plot_names = {"timeseries", "fft", "profile", "pulsestack", "tdfs", "modenv"};
  plot_dims  = [1, 1, 1, 2, 2, 1];
  plot_isfft = [0, 1, 0, 0, 1, 0];

  for n = 1:length(plot_names)
    figures.(plot_names{n}).drawfcn     = str2func(["plot_", plot_names{n}]); % Function that does the plotting
    figures.(plot_names{n}).fig_handle  = [];                                 % Handle to the figure
    figures.(plot_names{n}).ax_handle   = [];                                 % Handle(s) to the plot axes
    figures.(plot_names{n}).quantise    = 0;                                  % Whether a click is rounded to the nearest "image pixel"
    figures.(plot_names{n}).dims        = plot_dims(n);                       % Number of dimensions in plot
    figures.(plot_names{n}).isfft       = plot_isfft(n);                      % Whether the plot is of type "FFT"
  end % for

  % Set the default positions for the timeseries and fft figures
  screensize = get(0, 'screensize');
  gapsize    = 100;

  winsize_x  = screensize(3) - 2*gapsize;
  winsize_y  = floor((screensize(4) - 3*gapsize)/2);
  winpos_x   = gapsize;
  winpos_y   = 2*gapsize + winsize_y;
  figures.timeseries.defaultpos = [winpos_x, winpos_y, winsize_x, winsize_y];

  winpos_y   = gapsize;
  figures.fft.defaultpos = [winpos_x, winpos_y, winsize_x, winsize_y];

  % Load known pulsar periods
  load_periods();

  % Load initial file, if given
  if (nargin == 1)
    load_data([pwd(), "/"], init_filename);

    % Draw the main (=timeseries) window
    create_figure("timeseries");

    % Plot it up
    figures.timeseries.drawfcn();
  else
    % Just create the figure window
    create_figure("timeseries");
  end % if

  % Initially, ensure that there are no changes to be saved
  set_unsaved_changes(false);

end % function


function load_periods()

  global pulsars

  f = fopen('periods.dat');

  % First, count the lines
  nlines = fskipl(f, Inf);
  frewind(f);

  % Now read them in
  pulsars = struct();
  for n = 1:nlines
    [pulsarname, pulsarperiod] = fscanf(f, '%s %f', "C");
    pulsars.(pulsarname).period = pulsarperiod;
  end % for

end % function

%%%%%%%%%%%%%%%%%%%%%
% Everything after this point will eventually be relocated
%

function Xplot_modenv(src, data, fig_no)

  global fig_handles
  global plot_params

  global tdfs

  global period
  global only_visible

  % Switch to/Create 2DFS figure and keep track of the view window
  fig_no = 7;
  first_time = (fig_handles(fig_no) == 0);

  if (first_time)
    create_figure(0,0,fig_no);
  end % if

  figure(fig_handles(fig_no));

  % Assume: tdfs has already been calculated by now

  % Now calculate the inverse 2DFS (=m_{phi,t} in E&S(2002))
  M = ifft2(tdfs);

  % Do singular-value decomposition on m_{phi,t} = m_phi(phi) * m_t(t)
  [U,S,V] = svd(M);
  m_t     = U(:,1) * sqrt(S(1,1));
  m_phi   = V(:,1) * sqrt(S(1,1));

  % Calculate various x-axes
  if (only_visible)
    tmin = max([floor(plot_params(1,1)/period),0]) + 1;
    tmax = tmin + length(m_t) - 1;
  else
    tmin = 1;
    tmax = length(m_t);
  end % if

  pulse_no = [tmin:tmax];
  longitude = [0:(length(m_phi)-1)] / length(m_phi) * 360;

  subplot(2,2,1);
  plot(pulse_no, abs(m_t));
  xlabel("Pulse Number");
  ylabel("Amplitude");

  subplot(2,2,2);
  plot(longitude, abs(m_phi));
  xlabel("Pulse longitude");
  ylabel("Amplitude");

  subplot(2,2,3);
  plot(pulse_no, arg(m_t) * 180/pi, 'x');
  xlabel("Pulse Number");
  ylabel("Phase");

  subplot(2,2,4);
  plot(longitude, arg(m_phi) * 180/pi, 'x');
  xlabel("Pulse longitude");
  ylabel("Phase");

  % Calculate the "signal-to-noise" of the result
  SNR = S(1,1) / (trace(S) - S(1,1)); % <-- Not sure if this is the correct way to do it...

  set(fig_handles(fig_no), "Name", sprintf('Complex Envelope: \"SNR?\" = %f', SNR));

end % function

function Xclick_p2p3(src, data, fig_no)

  global fig_handles
  global quantise

  global orig_windowbuttondownfcn_p2p3
  global orig keypressfcn_p2p3

  h = fig_handles(fig_no);

  orig_windowbuttondownfcn_p2p3 = get(h, "windowbuttondownfcn");
  orig_keypressfcn_p2p3         = get(h, "keypressfcn");

  title_stem = ["Choose P2 and P3 (hat) by clicking on the 2DFS\nPress 'q' to turn quantisation"];
  title([title_stem, " ", on_off(~quantise)]);

  set(h, "windowbuttondownfcn", @collect_p2p3_click);
  set(h, "keypressfcn", {@collect_quantise_keypress, title_stem});

end % function

function Xcollect_p2p3_click(src, button)

  global period
  global P2hat
  global P3hat
  global timeseries_grid

  global orig_windowbuttondownfcn_p2p3
  global orig_keypressfcn_p2p3

  global quantise

  if (button == 1) % Left mouse click

    % Get the clicked point
    pos = get_curr_pos(gca, quantise);

    % Calculate P2 and P3 hat (remember to round for pixels)
    P2hat = period/pos(1);
    P3hat = period/pos(2);
    title("");
    set(gcf(), "windowbuttondownfcn", orig_windowbuttondownfcn_p2p3);
    set(gcf(), "keypressfcn", orig_keypressfcn_p2p3);

  end % if

  set_unsaved_changes(true);
  get_axes(0,0,6);
  replot(6);

end % function

