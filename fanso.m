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

function Xapply_shift_DC()

  global shift_DC
  global pos_relative
  global tdfs

  % First, a bit of error checking
  if (isempty(shift_DC))
    shift_DC = [0,0];
  end % if

  if (isempty(pos_relative))
    pos_relative = [0,0];
  end % if

  shift_amount = shift_DC + pos_relative;

  % Convert shift_amount to "pixel" units
  pixel_size = [1, 1/rows(tdfs)];
  shift_amount = round(shift_amount ./ pixel_size);

  % And do the shift (as simple as that!)
  tdfs = shift(tdfs, shift_amount(1), 2); % in the x-direction
  tdfs = shift(tdfs, shift_amount(2), 1); % in the y-direction

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

function Xcreate_filter(src, data, fig_no, hor_or_vert)
% hor_or_vert should be "horizontal" or "vertical"

  global fig_handles
  global new_filter
  global quantise

  global orig_windowbuttondownfcn_filter
  global orig_windowbuttonmotionfcn_filter
  global orig keypressfcn_filter

  hv = find(strcmp(hor_or_vert, {"horizontal", "vertical"})) - 1; % "h..." --> 0;  "v..." --> 1
  new_filter = [0, 0, hv, 0]; % [coord, width, hor_or_vel, marked_for_deletion]
  h = fig_handles(fig_no);

  orig_windowbuttondownfcn_filter   = get(h, "windowbuttondownfcn");
  orig_windowbuttonmotionfcn_filter = get(h, "windowbuttonmotionfcn");
  orig_keypressfcn_filter           = get(h, "keypressfcn");

  title_stem = ["Click the centre of the ", hor_or_vert, " filter\nPress 'q' to turn quantisation"];
  title([title_stem, " ", on_off(~quantise)]);

  set(h, "windowbuttondownfcn", {@collect_filter_click, 1, hor_or_vert}); % "1" = First click
  set(h, "windowbuttonmotionfcn", []); % <-- YET-TO-IMPLEMENT
  set(h, "keypressfcn", {@collect_quantise_keypress, title_stem});

  get_axes(0,0,fig_no);

end % function

function Xcollect_filter_click(src, button, click_no, hor_or_vert)

  global quantise
  global new_filter
  global filters

  global orig_windowbuttondownfcn_filter
  global orig_windowbuttonmotionfcn_filter
  global orig_keypressfcn_filter

  % Get the last clicked position
  pos = get_curr_pos(gca, quantise);
  idx = 2 - new_filter(3); % "h..." --> 2;  "v..." --> 1

  switch click_no
    case 1
      % Include the clicked value in the new filter
      new_filter(1) = pos(idx);

      % Change figure title
      title_stem = ["Click the edge of the ", hor_or_vert, " filter\nPress 'q' to turn quantisation"];
      title([title_stem, " ", on_off(~quantise)]);

      % Change callback functions on mouse click and key-press
      set(src, "windowbuttondownfcn", {@collect_filter_click, 2, hor_or_vert});
      set(src, "keypressfcn", {@collect_quantise_keypress, title_stem});
    case 2
      % Include the clicked value in the new filter and annex the new filter onto the list of filters
      new_filter(2) = abs(pos(idx) - new_filter(1));
      if (new_filter(2) > 0)
        no_filters = rows(filters);
        filters((no_filters+1),:) = new_filter;
        set_unsaved_changes(true);
      else
        errordlg("Cannot create filter of zero width");
      end % if

      % Reset new_filter to empty (so that nothing extraneous is plotted in plot_tdfs)
      new_filter = [];

      % Set things back to the way they were
      set(src, "windowbuttondownfcn", orig_windowbuttondownfcn_filter);
      set(src, "windowbuttonmotionfcn", orig_windowbuttonmotionfcn_filter);
      set(src, "keypressfcn", orig_keypressfcn_filter);
  end % switch

  replot(6);

end % function

function Xcollect_quantise_keypress(src, evt, newtitle_stem)

  if (evt.Character == 'q')

    global quantise
    quantise = ~quantise;

    title([newtitle_stem, " ", on_off(~quantise)]);

  end

end % function

function Xpos = get_curr_pos(h, quantised)

  point = get(h, "currentpoint");
  point = point(1,1:2);

  % If a figure handle was passed, then the "point" is in "figure"
  % coordinates (i.e. in pixels measured from lower-left corner of
  % plot window, and needs to be converted to axes coordinates
  if (isfigure(h))
    hax       = gca();
    ax_pos    = get(hax, "position");
    fig_pos   = get(h, "position");
    fig_width = fig_pos(3:4);
    ax_lims   = axis()(1:4);
    ax_width  = diff(reshape(ax_lims,2,2));
    ax_orig   = ax_lims(1,3);
    point     = (point./fig_width - ax_pos(1:2)) ./ ax_pos(3:4).*ax_width + ax_orig;
  end % if

  if (quantised)
    global timeseries_grid
    dx = 1;                               % Size of a pixel in the x-direction (in the 2DFS)
    dy = 1/rows(timeseries_grid);         % Size of a pixel in the y-direction (in the 2DFS)
    pos = [round(point(1)/dx) * dx, round(point(2)/dy) * dy];
  else
    pos = [point(1), point(2)];
  end % if

end % function

function Xdelete_filter(src, data, fig_no)

  global fig_handles
  global orig_windowbuttondownfcn_delfilter
  global orig_keypressfcn_delfilter
  global filters

  % "Save" the old
  h = fig_handles(fig_no);
  orig_windowbuttondownfcn_delfilter = get(h, "windowbuttondownfcn");

  set(h, "windowbuttondownfcn", @collect_deletefilter_click);
  set(h, "keypressfcn", @keypress_deletefilter);

  title("Left click = select/deselect\nRight click = cancel\n'd' = delete")

  get_axes(0,0,fig_no);

end % function

function Xcollect_deletefilter_click(src, button)

  global filters
  global shift_DC
  global orig_windowbuttondownfcn_delfilter
  global orig_keypressfcn_delfilter

  switch button
    case 1 % Left click = mark for deletion
      % Get current position of cursor and adjust by shift_DC amount
      pos = get_curr_pos(gca, false) - shift_DC;

      % Select filters containing cursor
      h_distances = (filters(:,3) == 0) .* (abs(pos(2) - filters(:,1)) ./ filters(:,2)); % horizontal filters
      v_distances = (filters(:,3) == 1) .* (abs(pos(1) - filters(:,1)) ./ filters(:,2)); % vertical filters
      distances = h_distances + v_distances;
      [mindist, closest] = min(distances);

      if (mindist <= 1)
        filters(closest,4) = ~filters(closest,4);
        replot(6);
        title("Left click = select/deselect\nRight click = cancel\n'd' = delete")
      end % if

    case 3 % Right click = cancel
      filters(:,4) = 0;
      if (isempty(filters))
        clear filters
      end % if

      set(src, "windowbuttondownfcn", orig_windowbuttondownfcn_delfilter);
      set(src, "keypressfcn", orig_keypressfcn_delfilter);
      title("");
      replot(6);
  end % switch

end % function

function Xkeypress_deletefilter(src, evt)

  global filters
  global orig_windowbuttondownfcn_delfilter
  global orig_keypressfcn_delfilter

  if (evt.Character == 'd')
    to_be_deleted = filters(:,4);
    if (any(to_be_deleted))
      set_unsaved_changes(true);
      filters = filters(~filters(:,4),:);
    end % if

    if (isempty(filters))
      clear filters
    end % if

    set(src, "windowbuttondownfcn", orig_windowbuttondownfcn_delfilter);
    set(src, "keypressfcn", orig_keypressfcn_delfilter);
    replot(6);
  end % if

end % function

function Xclick_shiftDC(src, data, fig_no)

  global fig_handles
  global orig_windowbuttondownfcn_shiftDC

  h = fig_handles(fig_no);
  orig_windowbuttondownfcn_shiftDC = get(h, "windowbuttondownfcn");

  get_axes(fig_no);

  set(h, "windowbuttondownfcn", @buttondown_shiftDC);

  title("Click and drag to shift the origin\nRight click to confirm");

end % function

function Xbuttondown_shiftDC(src, button)

  switch button
    case 1 % Left button = initiate drag
      global pos_start
      global orig_windowbuttonmotionfcn_shiftDC
      global orig_windowbuttonupfcn_shiftDC

      pos_start = get_curr_pos(gcf, true);

      orig_windowbuttonmotionfcn_shiftDC = get(src, "windowbuttonmotionfcn");
      orig_windowbuttonupfcn_shiftDC     = get(src, "windowbuttonupfcn");

      set(src, "windowbuttonmotionfcn", @buttonmotion_shiftDC)
      set(src, "windowbuttonupfcn",     @buttonup_shiftDC);
    case 3 % Right button = confirm and stop
      global orig_windowbuttondownfcn_shiftDC
      set(src, "windowbuttondownfcn", orig_windowbuttondownfcn_shiftDC);
      title("");
  end % switch
  
end % function

function Xbuttonmotion_shiftDC(src, button)

  global pos_start
  global pos_relative % <-- this global variable will get used in plot_tdfs()
                      %     for "live" shifting feedback

  % Get current (quantised) mouse position
  pos_curr = get_curr_pos(gcf, true);

  % Work out how much has been moved since the initial button down
  pos_relative = pos_curr - pos_start;

  % Update plot
  replot(6);
  title("Click and drag to shift the origin\nRight click to confirm");

end % function

function Xbuttonup_shiftDC(src, button)

  global pos_relative
  global shift_DC

  global orig_windowbuttonmotionfcn_shiftDC
  global orig_windowbuttonupfcn_shiftDC

  % Update shift_DC and reset pos_relative
  if (any(pos_relative ~= 0))
    set_unsaved_changes(true);
    shift_DC        += pos_relative;
    pos_relative(:)  = 0;
  end % if

  set(src, "windowbuttonmotionfcn", orig_windowbuttonmotionfcn_shiftDC);
  set(src, "windowbuttonupfcn",     orig_windowbuttonupfcn_shiftDC);

end % function

