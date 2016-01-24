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
  disp('This function requires ',gtk,' to be installed.')
  return;
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

m = uimenu("label", "FF&T", "accelerator", "F", "callback", {@plot_fft, dt, timeseries});
%plot_fft(timeseries)

end % function


function plot_fft(src,data,dt,timeseries)
% This is a callback function
  % Set position,size of figure 2 window
  global fig_x
  global fig_y
  global defaultwidth
  global defaultheight
  figure(2, "Position", [fig_x(2), fig_y(2), defaultwidth, defaultheight]);

  % Get only portion of timeseries visible in Figure 1
  figure(1);
  a = axis();
  min_idx = max([floor(a(1)/dt)+1, 1]);
  max_idx = min([floor(a(2)/dt)+1, length(timeseries)]); % <-- Check this
  timeseries = timeseries([min_idx:max_idx]);

  % Plot it up
  figure(2);
  ffted = fft(timeseries);
  absed = abs(ffted);
  plot(absed);
end % function


