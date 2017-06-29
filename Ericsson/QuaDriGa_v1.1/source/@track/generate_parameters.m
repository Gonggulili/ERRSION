function par = generate_parameters( h_track, overlap, usage, check_parfiles, verbose )
%GENERATE_PARAMETERS Generate large scale parameters and store them in "par"
% 
% GENERATE_PARAMETERS extracts the large scale parameters (LSPs) for the
% given scenario from the "parameter_set" class and stores them in
% "track.par". Hence, it automatically generates the LSPs and, thus,
% implements an easy-to-use interface for the "parameter_set" class.
%
% Since the track class does not handle transmitter positions, a default
% position of [0,0,25] is assumed. Please refer to
% "layout.generate_parameters" for a more detailed description.
%
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


% Parse Input variables
if exist( 'verbose' , 'var' )
    if ~( isnumeric(verbose) || ~isnan(verbose) || all(size(verbose) == [1 1]) )
        error('"verbose" must be numeric.')
    else
        verbose = logical(verbose);
    end
else
    verbose = true;
end

if ~exist( 'overlap' , 'var' )
    overlap = 0.5;
end

if ~exist( 'usage' , 'var' )
    usage = 2;
end

if ~exist( 'check_parfiles' , 'var' )
    check_parfiles = 1;
end

% Create a layout object that contains the track
h_layout = layout;
h_layout.simpar.show_progress_bars = verbose;
h_layout.track = h_track;
h_layout.tx_position(3) = 25;

% Call the parameter generation procedure for the layout object
tmp = h_layout.generate_parameters( overlap, usage , check_parfiles );
par = tmp{1};

end
