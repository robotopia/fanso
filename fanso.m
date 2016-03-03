% -- Program: fanso
%     Analyse a timeseries using the Fourier ANalysis Suite for Octave.
%     Requires FLTK.
%
%     Takes one optional argument, a .fan "FANSO" file.
%
%     Written by Sam McSweeney, 2016, Creative Commons Licence
%     sammy.mcsweeney@gmail.com

% First, check if there is already an instance of FANSO open
global __FANSO_INSTANCE__
if (isempty(__FANSO_INSTANCE__)) % i.e. there is NO instance already open
  __FANSO_INSTANCE__ = 1;
else
  % Exit gracefully
  disp("FANSO already running");
  return
end

% Second, check that Octave 4.0.0 is being used
v = version();
if (str2num(v(1)) < 4)
  disp("This software requires at least 4.0.0");
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

plot_names  = {"timeseries", "fft", "profile", "pulsestack", "tdfs", "modenv"};
plot_dims   = [1, 1, 1, 2, 2, 1];
plot_isfft  = [0, 1, 0, 0, 1, 0];
plot_nplots = [1, 1, 1, 1, 1, 4];

for n = 1:length(plot_names)
  figures.(plot_names{n}).drawfcn     = str2func(["plot_", plot_names{n}]); % Function that does the plotting
  figures.(plot_names{n}).fig_handle  = [];                                 % Handle to the figure
  figures.(plot_names{n}).ax_handle   = [];                                 % Handle(s) to the plot axes
  figures.(plot_names{n}).quantise    = 0;                                  % Whether a click is rounded to the nearest "image pixel"
  figures.(plot_names{n}).dims        = plot_dims(n);                       % Number of dimensions in plot
  figures.(plot_names{n}).isfft       = plot_isfft(n);                      % Whether the plot is of type "FFT"
  figures.(plot_names{n}).nplots      = plot_nplots(n);                     % Number of subplots in this figure
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

% Load initial file, if given
if (nargin >= 1)
  arg_list = argv();
  init_filename = arg_list{1};
  load_data("", init_filename); % <--- cause for errors here, if user supplies absolute path

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


