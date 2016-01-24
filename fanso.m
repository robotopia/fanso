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
global global_dt         = dt
global global_timeseries = timeseries

% Calculate default dimensions of figure windows
global screensize
global gapsize
global winsize_x
global winsize_y
global winpos_x
global winpos_y
global fig1 = figure();
global fig2 = figure();

global zeromean;      % 0 = Do nothing;               1 = Zero mean before applying FFT
global only_visible;  % 0 = FFT of entire timeseries; 1 = FFT of only visible timeseries
global apply_hamming; % 0 = Do nothing;               1 = Apply Hamming window
global apply_hanning; % 0 = Do nothing;               1 = Apply Hanning window

% Set initial values
screensize = get(0, 'screensize');
gapsize    = 100;

zeromean      = 0;
only_visible  = 0;
apply_hamming = 0;
apply_hanning = 0;

% Require that we are using the correct graphics toolkit
gtk = 'fltk';
if (~any(strcmp(available_graphics_toolkits(), gtk)))
  disp('This function requires ',gtk,' to be installed.')
  return;
else
  graphics_toolkit(gtk);
end % if

% Set up the menus
figure(fig1);
m_fft                = uimenu('label', 'FF&T', 'accelerator', 't');
m_fft_window         = uimenu(m_fft, 'label', '&Windowing function', 'accelerator', 'w');
m_fft_window_hamming = uimenu(m_fft_window, 'label', 'Ha&mming', 'accelerator', 'm', 'callback', @toggle_hamming);
m_fft_window_hanning = uimenu(m_fft_window, 'label', 'Ha&nning', 'accelerator', 'n', 'callback', @toggle_hanning);
m_fft_zeromean       = uimenu(m_fft, 'label', '&Zero-mean', 'accelerator', 'z', 'callback', @toggle_zeromean);
m_fft_visible        = uimenu(m_fft, 'label', 'Only &visible', 'accelerator', 'v', 'callback', @toggle_visible);

% Replot FFT when timeseries is zoomed
addlistener(gca, 'xlim', @plot_fft);

% Draw the plots for the first time
plot_timeseries(0, 0, 1);
plot_fft(0, 0, 1);

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

function plot_timeseries(src, data, first_time = 0)

  global fig1
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

  global global_timeseries
  global global_dt

  % Switch to timeseries figure and keep track of the view window
  if (first_time)
    figure(fig1, "Position", [winpos_x, winpos_y, winsize_x, winsize_y]);
  else
    figure(fig1);
    ax = axis();
  end % if

  % Calculate the values for the timeseries abscissa
  N = length(global_timeseries);
  t = [0:(N-1)] * global_dt;

  % Plot the timeseries!
  figure(fig1);
  plot(t, global_timeseries);
  xlabel('Time (s)');
  ylabel('Timeseries values');
  if (~first_time)
    axis(ax);
  end % if

end % function

function plot_fft(src, data, first_time = 0)
  % Set position,size of figures window
  global fig1
  global fig2

  global screensize
  global gapsize
  global winsize_x
  global winsize_y
  global winpos_x
  global winpos_y

  winsize_x  = screensize(3) - 2*gapsize;
  winsize_y  = floor((screensize(4) - 3*gapsize)/2);
  winpos_x   = [gapsize, gapsize];
  winpos_y   = [2*gapsize + winsize_y, gapsize];

  global zeromean
  global only_visible
  global apply_hamming
  global apply_hanning

  global global_timeseries
  global global_dt

  % Switch to FFT figure and keep track of the view window of timeseries figure
  figure(fig1, "Position", [winpos_x(1), winpos_y(1), winsize_x, winsize_y]);
  ax1 = axis();
  if (first_time)
    figure(fig2);
  else
    figure(fig2, "Position", [winpos_x(2), winpos_y(2), winsize_x, winsize_y]);
    ax2 = axis();
  end % if

  % Are we processing the whole timeseries or just the visible part?
  if (only_visible)
    min_idx = max([floor(ax1(1)/global_dt)+1, 1]);
    max_idx = min([floor(ax1(2)/global_dt)+1, length(global_timeseries)]); % <-- Check this
    to_be_ffted = global_timeseries([min_idx:max_idx]);
  else
    to_be_ffted = global_timeseries;
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
  plot(f, power);
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


