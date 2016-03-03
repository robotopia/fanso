function par_vals = read_par(par_filename)
%     Read in a ephemeris file and retrieve the required quanitites.
%     (Currently assumes a PSRCAT-like format)


% Read in the emphemeris into a cell array
par_data = importdata(par_filename, delimiter='\t');

% loop through each "row" of the ephemeris and assign the fields to a scalar structure
for n = 1:length(par_data)
  row = strsplit(strtrim(par_data{n}));
  par_vals.(row{1}) = row{2};


end % function
