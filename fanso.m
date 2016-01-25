function fanso(timeseries, dt)
% -- Function: fanso (timeseries, dt)
%     Analyse a timeseries using the Fourier ANalysis Suite for Octave.
%     Requires FLTK.
%
%     timeseries = a column vector of real values.
%
%     dt = a scalar (representing the time between consecutive
%          points in timeseries).
%
%     Written by Sam McSweeney, 2016, Creative Commons Licence
%     sammy.mcsweeney@gmail.com

% Create global variables which will be used for the actual plotting
global global_dt         = dt;
global global_timeseries = timeseries;
global global_flattened  = timeseries; % timeseries flattened by linear detrending between breakpoints
global breakpoint_mask   = ones(size(timeseries)); % 0 = do not use this value; 1 = use this value

% Variables for the size and position of the figure windows
global screensize;
global gapsize;
global winsize_x;
global winsize_y;
global winpos_x;
global winpos_y;
global fig1;
global fig2;

% Variables for the plot settings
global zeromean;      % 0 = Do nothing;               1 = Zero mean before applying FFT
global only_visible;  % 0 = FFT of entire timeseries; 1 = FFT of only visible timeseries
global apply_hamming; % 0 = Do nothing;               1 = Apply Hamming window
global apply_hanning; % 0 = Do nothing;               1 = Apply Hanning window
global apply_bps;     % 0 = Do nothing;               1 = Apply breakpoints (i.e. "flatten" timeseries)

% A global variable for the breakpoints
global breakpoints;

% Set initial values
screensize = get(0, 'screensize');
gapsize    = 100;

zeromean      = 0;
only_visible  = 0;
apply_hamming = 0;
apply_hanning = 0;
apply_bps     = 0;

breakpoints = [];

% Require that we are using the correct graphics toolkit
gtk = 'fltk';
if (~any(strcmp(available_graphics_toolkits(), gtk)))
  disp('This function requires ',gtk,' to be installed.')
  return;
else
  graphics_toolkit(gtk);
end % if

% Draw the plots for the first time
plot_timeseries();

end % function




function toggle_hamming(src, data)

  global apply_hamming
  apply_hamming = ~apply_hamming;

  if (apply_hamming)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  plot_fft()

end % function

function toggle_hanning(src, data)

  global apply_hanning
  apply_hanning = ~apply_hanning;

  if (apply_hanning)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  plot_fft()

end % function

function toggle_zeromean(src, data)

  global zeromean
  zeromean = ~zeromean;

  if (zeromean)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  plot_fft()

end % function

function toggle_visible(src, data)

  global only_visible
  only_visible = ~only_visible;

  if (only_visible)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  plot_fft()

end % function

function add_breakpoints(src, data)

  global fig1;
  global breakpoints;
  global apply_bps;
  global m_bp_apply;

  if (apply_bps)
    toggle_flatten(m_bp_apply);
  end % if

  do
    figure(fig1);
    title("Left mouse button = add breakpoint; right = remove breakpoint; 's' = stop\nDON'T CLOSE THIS WINDOW!");
    [x, y, button] = ginput(1);
    switch button
      case 1 % left mouse button
        breakpoints = union(breakpoints, x);
      case 3 % right mouse button
        if (length(breakpoints) > 0)
          [nearest_bpdiff, nearest_idx] = min(abs(breakpoints - x));
          breakpoints = setdiff(breakpoints, breakpoints(nearest_idx));
        end
    end % switch
    flatten();
    plot_timeseries(0,0,0);
  until (button == 115) % 115='s'
  title('');

end % function

function toggle_flatten(src, data)

  global fig2

  global apply_bps
  apply_bps = ~apply_bps;

  if (apply_bps)
    set(src, 'checked', 'on');
  else
    set(src, 'checked', 'off');
  end % if

  % Redraw plots
  plot_timeseries();

  if (~isempty(fig2))
    if (isfigure(fig2))
      plot_fft();
    end % if
  end % if

end % function

function plot_timeseries(src, data)

  global fig1

  global global_timeseries
  global global_flattened
  global global_dt

  global breakpoints
  global apply_bps

  % Switch to/Create timeseries figure and keep track of the view window
  first_time = 0;

  if (isempty(fig1))
    first_time = 1;
  end

  if (~isfigure(fig1))
    first_time = 1;
  end

  if (first_time)

    global screensize
    global gapsize
    global winsize_x
    global winsize_y
    global winpos_x
    global winpos_y

    winsize_x  = screensize(3) - 2*gapsize;
    winsize_y  = floor((screensize(4) - 3*gapsize)/2);
    winpos_x   = gapsize;
    winpos_y   = 2*gapsize + winsize_y;

    fig1 = figure("Position", [winpos_x, winpos_y, winsize_x, winsize_y]);

    % Set up the menus
    m_fft                = uimenu('label', 'FF&T');
    m_fft_window         = uimenu(m_fft, 'label', '&Windowing function', 'accelerator', 'w');
    m_fft_window_hamming = uimenu(m_fft_window, 'label', 'Ha&mming', 'accelerator', 'm', 'callback', @toggle_hamming);
    m_fft_window_hanning = uimenu(m_fft_window, 'label', 'Ha&nning', 'accelerator', 'n', 'callback', @toggle_hanning);
    m_fft_zeromean       = uimenu(m_fft, 'label', '&Zero-mean', 'accelerator', 'z', 'callback', @toggle_zeromean);
    m_fft_visible        = uimenu(m_fft, 'label', 'Only &visible', 'accelerator', 'v', 'callback', @toggle_visible);
    m_fft_plotfft        = uimenu(m_fft, 'label', 'Plot FFT', 'separator', 'on', 'callback', @plot_fft);
    m_bp                 = uimenu('label', '&Breakpoints');
    m_bp_add             = uimenu(m_bp, 'label', 'Add &breakpoints', 'accelerator', 'b', 'callback', @add_breakpoints);
    global m_bp_apply    = uimenu(m_bp, 'label', 'A&pply breakpoints', 'separator', 'on', 'callback', @toggle_flatten);

  else
    figure(fig1);
    ax = axis();
  end % if

  % Calculate the values for the timeseries abscissa
  N = length(global_timeseries);
  t = [0:(N-1)]' * global_dt;

  % Plot different things depending on whether the breakpoints are
  % to be "applied" or not (i.e. whether the timeseries has been flattened).
  if (apply_bps)
    plot(t, global_flattened, 'b');
  else
    % Prepare vertical lines for breakpoint plotting
    if(~isempty(breakpoints))
      bp_xs = [1;1] * breakpoints;
      ymin = min(global_timeseries);
      ymax = max(global_timeseries);
      bp_ys = repmat([ymin;ymax],size(breakpoints));

      % Plot the flattened timeseries, the original timeseries, and breakpoints
      plot(t, global_flattened,  'g', ...
           t, global_timeseries, 'b', ...
           bp_xs, bp_ys, 'r', 'linewidth', 2.0);
    else
      % Plot the timeseries!
      plot(t, global_timeseries, 'b');
    end % if
  end % if
  figure(fig1)
  xlabel('Time (s)');
  ylabel('Timeseries values');
  if (~first_time)
    axis(ax);
  end % if

end % function

function plot_fft(src, data)
  % Set position,size of figures window
  global fig1
  global fig2

  global global_timeseries
  global global_flattened
  global global_dt

  global zeromean
  global only_visible
  global apply_hamming
  global apply_hanning
  global apply_bps

  % Make sure there is a timeseries figure showing
  if (isempty(fig1))
    error("No timeseries plot found");
  end % if

  if (~isfigure(fig1))
    error("No timeseries plot found");
  end % if

  figure(fig1);
  ax1 = axis();

  % Set up the FFT figure
  first_time = 0;
  if (isempty(fig2))
    first_time = 1;
  end

  if (~isfigure(fig2))
    first_time = 1;
  end

  if (first_time)

    global screensize
    global gapsize
    global winsize_x
    global winsize_y
    global winpos_x
    global winpos_y

    winsize_x  = screensize(3) - 2*gapsize;
    winsize_y  = floor((screensize(4) - 3*gapsize)/2);
    winpos_x   = gapsize;
    winpos_y   = gapsize;

    fig2 = figure("Position", [winpos_x, winpos_y, winsize_x, winsize_y]);
  else
    figure(fig2);
  end % if
  ax2 = axis();

  % Are we plotting the original timeseries or the flattened timeseries?
  if (apply_bps)
    to_be_ffted = global_flattened;
  else
    to_be_ffted = global_timeseries;
  end % if

  % Are we processing just the visible part?
  if (only_visible)
    min_idx = max([floor(ax1(1)/global_dt)+1, 1]);
    max_idx = min([floor(ax1(2)/global_dt)+1, length(global_timeseries)]); % <-- Check this
    to_be_ffted = to_be_ffted([min_idx:max_idx]);
  end % if

  % Get the length of the timeseries to be FFT'd
  n  = length(to_be_ffted);

  % Apply zeromean
  if (zeromean)
    to_be_ffted = to_be_ffted - mean(to_be_ffted);
  end % if

  % Apply Hamming window
  if (apply_hamming)
    to_be_ffted = hamming(n) .* to_be_ffted;
  end % if

  % Apply Hanning window
  if (apply_hanning)
    to_be_ffted = hanning(n) .* to_be_ffted;
  end % if

  % Calculate the values...
  % ...for the FFT abscissa...
  df = 1/(global_dt*n);
  f  = [0:(n-1)] * df;
  % ...and the FFT ordinate (=power)...
  ffted = fft(to_be_ffted);
  absed = abs(ffted);
  power = absed.^2;

  % Plot up the FFT
  figure(fig2);
  plot(f, power, 'b');
  xlabel('Frequency (Hz)');
  ylabel('Power');
  if (~first_time)
    % Keep xlim's as they were, but rescale ylim's accordingly
    idx = find((f >= ax2(1)) & (f <= ax2(2)));
    ymax = max(power(idx));
    ax2 = [ax2(1), ax2(2), 0, ymax];
    axis(ax2);
  end % if

end % function


function flatten()

  global global_dt
  global global_timeseries
  global global_flattened
  global breakpoint_mask

  global breakpoints

  if (isempty(breakpoints))
    global_flattened = global_timeseries;
  end % if

  % Book-end the breakpoints with initial and final values
  N = length(global_timeseries);
  t = [0:(N-1)]' * global_dt;

  bps = [t(1)-1, breakpoints, t(end)+1];

  % Loop through the breakpoints and detrend each segment
  for i = 1:(length(bps)-1)
    flatten_idxs  = (t >= bps(i)) & (t < bps(i+1));
    model_idxs    = flatten_idxs & breakpoint_mask;
    fit           = polyfit(t(model_idxs), global_timeseries(model_idxs), 1);
    m             = fit(1);
    c             = fit(2);
    y_orig        = global_timeseries(flatten_idxs);
    x             = t(flatten_idxs);
    global_flattened(flatten_idxs) ...
                  = y_orig - (m*x + c);
  end % for

end % function
