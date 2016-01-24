function fanso(timeseries, dt)

% Calculate default dimensions of figure windows
global screensize = get(0,'screensize');
global gapsize = 100; % pixels
global defaultwidth = screensize(3) - 2*gapsize;
global defaultheight = floor((screensize(4) - 3*gapsize)/2);
global fig_x = [gapsize, gapsize];
global fig_y = [2*gapsize + defaultheight, gapsize];


% Require that we are using the correct graphics toolkit
gtk = 'fltk';
if (~any(strcmp(available_graphics_toolkits(), gtk)))
  disp('This program requires ',gtk,' to be installed.')
  exit
else
  graphics_toolkit(gtk);
end % if

% Get max and min values in timeseries
ymax = max(timeseries);
ymin = min(timeseries);

% Plot the timeseries!
figure(1, "Position", [fig_x(1), fig_y(1), defaultwidth, defaultheight]);
N = length(timeseries);
t = [0:(N-1)] * dt;

plot(t, timeseries);
xlabel('Time (s)');
ylabel('Timeseries values');

m = uimenu("label", "FF&T", "accelerator", "F", "callback", {@plot_fft, "timeseries"});

end % function


% Usage function
function my_usage()
  disp('usage: fanso FILENAME')
  disp('  FILENAME must be an Octave-produced ASCII file containing')
  disp('  the variables ''dt'', a scalar representing the time between')
  disp('  consecutive data points (in seconds), and ''timeseries'', the')
  disp('  values of the timeseries data points.')
  exit
end % function

% Display online help -- Prefer to user menus and accelerators. Delete me when this function is no longer used.
function online_help()
  disp('b\tToggle breakpoints mode (init = off). When on, left mouse button')
  disp('\tclicks in the timeseries window will introduce a breakpoint. Right')
  disp('\tclicks will remove the nearest breakpoint.')
  disp('h\tShows this online help.')
end % function

function plot_fft(timeseries)
timeseries
  global fig_x
  global fig_y
  global defaultwidth
  global defaultheight
  figure(2, "Position", [fig_x(2), fig_y(2), defaultwidth, defaultheight]);
  ffted = fft(timeseries);
  absed = abs(ffted);
  plot(absed);
end % function


