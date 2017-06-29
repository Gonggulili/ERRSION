%% Generation of Satellite Channels
%
% This tutorial introduces a simple automatic state sequence generator
% (SSG), it shows how to set satellite orbital positions, adjust parameters
% and generate a time series of channel coefficients.


%% Setting up the Simulation Parameters
% First, we set up the general simulation parameters. We choose a center
% frequency of 2.1 GHz. We also want to use drifting in order to get the
% correct delays and angles for the time-continuous simulation. A sample
% density of 2.5 ensures that the channel coefficients can be interpolated
% to different playback speeds later on.

close all
clear all

set(0,'defaultTextFontSize', 14)    % Set default font size for the plots
set(0,'defaultAxesFontSize', 14)

s = simulation_parameters;          % Basic simulation parameters
s.center_frequency      = 2.185e9;
s.samples_per_meter     = 2;


%% Creating a random User Trajectory
% Next, we generate a simulation track. A track describes the  movement of
% a mobile terminal. It is composed of an ordered list of positions. During
% the simulation, one snapshot is generated for each position on the track.
% Later on, the generation of the track is done by the state sequence
% generator. Here, we implement a simple version of the sequence generator
% to generate a random track.

%%
% We first create a set of streets with different length. We assume a
% normal distribution of the street length where the parameters mu and
% sigma were fitted from random distances between two crossings in central
% Berlin (measured with Google earth).

street_length_mu = 187;         % Average street length in m
street_length_sigma = 83;
min_street_length = 50;

turn_probability = 0.5;         % The prob. that the car turns at a crossing
curve_radius = 10;              % The curve radius in m

%%
% For the given parameters, we calculate a list of points along the
% track that resemble the street grid and the turns at crossings.

diro  = 0;                      % Initial start diretion
point = 0;                      % The start point (always at [0,0])
m = 1;                          % A counter for the points
for n = 1:12                    % We simulate 20 street segments

    % Get a random street length drawn from the distribution defined above
    street_length = randn*street_length_sigma + street_length_mu;
    while street_length < min_street_length
        street_length = randn*street_length_sigma + street_length_mu;
    end
    
    % Get 3 points along the street
    point(m+1) = point(m) + exp(1j*diro) * street_length*0.1;
    point(m+2) = point(m) + exp(1j*diro) * street_length*0.9;
    point(m+3) = point(m) + exp(1j*diro) * street_length;
    m=m+3;
    
    % At a crossing, the car could change its direction. This is
    % modeled here
    if rand < turn_probability
        
        dirb = pi;
        while dirb >= pi/2
            dirn = diro + sign( rand-0.5 ) * pi/2 + randn*pi/12;
            dirb = abs( angle( exp(1j*dirn) ));
        end
        
        point(m+1) = point(m) + curve_radius*( exp(1j*diro) + exp(1j*dirn) );
        diro = dirn;
        m=m+1;
    end
end
point = point .* exp(1j*rand*2*pi);     % Random start direction

%%
% Next, we create a track object and pass the points along the track. We then 
% use the internal interpolation functions to interpolate the track to 1 
% point per meter.

t = track;                      % Create a track object
t.positions = [ real(point) ; imag(point) ; zeros(1,numel(point))];
t.interpolate_positions( 1 );   % Interpolate to 1 point per meter

t.visualize;                    % Plot Track


%% Define state sequence and segments
% We now assemble a rudimentary state sequence generator that generates
% different states along the track. We first define the distribution
% parameters of the segment length and then calculate the segments
% themselves. The two possible states are "MIMOSA_10-45_LOS" which stands 
% for LOS or good state and "MIMOSA_10-45_NLOS" for NLOS or bad state.

segment_length_mu = 30;         % Average segment length in m
segment_length_sigma = 12;      % Standard deviation in m
min_segment_length = 10;        % Minimum segment length in m

% Now we define the segments (the states) along the track
ind = 1;
while ind < t.no_snapshots

    % Each scenario has a 50% probability
    if rand < 0.5
        t.scenario{ t.no_segments } = 'MIMOSA_10-45_LOS' ;
    else
        t.scenario{ t.no_segments } = 'MIMOSA_10-45_NLOS' ;
    end
    
    % Get the length of the current segment
    segment_length = randn*segment_length_sigma + segment_length_mu;
    while segment_length<min_segment_length
        segment_length = randn*segment_length_sigma + segment_length_mu;
    end
    segment_length = round(segment_length);     % Segment length
    ind = ind + segment_length;                 % Start of next segment
    
    if ind < t.no_snapshots     % Exception for the last segment
        t.no_segments = t.no_segments + 1;
        t.segment_index( t.no_segments ) = ind;
    end
end

%%
% Finally, we interpolate the track to the given sample density (2 samples
% per half-wave-length) and plot the track.

t.interpolate_positions( s.samples_per_meter );
t.visualize;


%% Defining Antenna Arrays
% In the third step, we set up our antenna arrays for the transmitter at
% the satellite and the receiver. We use  synthetic dipole antennas for
% this case. Two dipoles are crossed by an angle of 90 degree. The signal is then
% split and fed with a 90 degree phase shift to both elements generating RHCP and
% LHCP signals.

% Create a patch antenna with 120 degree opening
a = array('custom',120,120,0);

% Copy element 1 to element 2 - the resulting antenna array has two
% elements, both dipoles.
a.copy_element(1,2);

% Rotate the second pattern by 90 degree around the x-axis.
a.rotate_pattern(90,'x',2);

%%
% Set the coupling between the elements. The Tx-signal for the first
% element is shifted by +90 degree out of phase and put on the second element.
% The signal for the second element is shifted by -90 degree and copied to the
% first element. Both antennas thus radiate a RHCP and a LHCP wave.

a.coupling = 1/sqrt(2) * [1 1;1j -1j];      % LHCP and RHCP

% Create a copy of the array for the receiver.
b = a.copy_objects;

% Rotate the receive antenna array to face sky-wards.
b.rotate_pattern(-90,'y');

b.visualize;                    % Plot the pattern of the Rx-Antenna


%% Setting up the Layout
% In this step, we combine the track, the antennas and the position of the
% satellite into a simulation layout. A layout object contains all the
% geometric information that are necessary to run the simulation. First,
% we define the position of the satellite. Since the model uses Cartesian
% coordinates, we have to transform the position of the satellite first.

l = layout( s );                % Create a new layout

% Choose a random satellite position (Astra 2, seen from Berlin).
% The distance only needs to be big enough to ensure insignificant changes
% in the reception angle on the ground.

sat_el      = 31.6;             % Elevation angle
sat_az      = 180;              % Azimuth angle (South = 180 degree)
rx_latitude = 53.5;             % Latitude of the Rx

% Approximate the satelite distance for GEO orbit
dist_x      = 35786 + rx_latitude/90 * 6384;    % [km]
dist_y      = (1-rx_latitude/90) * 6384;        % [km]
sat_dist    = sqrt(dist_x^2 + dist_y^2);        % [km]
sat_dist    = sat_dist*1e3;                     % [m]

% Transform angles to Cartesian coordinates
sat_x = sat_dist * cosd(sat_el) * cosd( -sat_az+90 );
sat_y = sat_dist * cosd(sat_el) * sind( -sat_az+90 );
sat_z = sat_dist * sind(sat_el);

% We also turn the antenna of the satellite so it points to the receiver.
a.rotate_pattern( sat_el , 'y' );
a.rotate_pattern( 270-sat_az , 'z' );

% Set the satellite position in the layout
l.tx_position = [ sat_x ; sat_y ; sat_z ];

l.track = t;                    % Set the track for the receiver
l.tx_array = a;                 % Set the tx_array
l.rx_array = b;                 % Set the rx_array


%% Setting up scenario parameters
% Next, the large scale parameters are set. The first line calls
% "l.create_parameter_sets", a built-in function that processes the data in
% the layout and returns a new "parameter_set" object "p". "p" is an array
% with two elements. One of them contains all the parameters for the good
% state (LOS) and one for the bad state (NLOS).  

p = l.create_parameter_sets(0);

%%
% Each parameter set has two different kinds of parameters. One for the
% scenario and one for the current state. For example, a scenario might
% have an average RMS Delay spread of 158 ns plus a certain variance which
% defines a range for the RMSDS. In addition to that, there are
% cross-correlations with other parameters such as the angular spread at
% the transmitter. All those parameters are stored in the "scenpar"
% property. For the good state, that parameters are:

Sl = strcmp( { p(1).name , p(2).name } ,'MIMOSA-10-45-LOS_Tx01' );  % Select good state
p(Sl).scenpar                                        % Show parameter list

%%
% Note here, that the values are given for a log-normal distribution. Thus,
% the RMSDS in nanoseconds follows from

10^( p(Sl).scenpar.DS_mu ) * 1e9

%%
% Each parameter on that list can be changed by just assigning it a new
% value. Here, we set the number of clusters for the LOS scenario to 7.
% Note that the default settings are stored in files in the sub-folder 
% "config" of the channel model folder. Here, the default settings can be 
% permanently set.
%
% After a change, the parameters of the segments need to be updated. This
% is done by calling the "update_parameters" method.

p(Sl).scenpar.NumClusters = 7;
p.update_parameters;

%%
% When "update_parameter" is called, the specific parameters for each
% segment are generated. E.g. each segment gets assigned a RMS Delay Spread
% and other values which are drawn from the statistics defined in scenpar.
% For the LOS segments, the individual RMSDS values for each segment are:

rmsds = p(Sl).ds*1e9
average = mean(p(Sl).ds*1e9)


%% Generate channel coefficients
% Next, we generate the channel coefficients. This is a lengthy task. The
% next line then combines the channels of the individual segments into a
% time-continuous channel. Here, the parameter (0.2) sets the length of the
% overlap region between two segments. In this case, it is 20%. 

c = p.get_channels;                             % Generate coefficients
cn = c.merge(0.2);                              % Combine segments


%% Evaluation of the data
% The next two plots show some basic evaluations of the generated
% coefficients. The first plot shows the received power for the 4 MIMO
% links along the track. The plot shows the differences between the LOS and
% NLOS segments and the cross-pol discrimination between the MIMO links.
% The average path loss for LOS was set to -95 dB and for NLOS -113 dB.

dist = (1:cn.no_snap)*t.get_length/cn.no_snap;
ind  = find(strcmp(t.scenario,'MIMOSA_10-45_LOS'));
los  = [];
for n = 1:numel(ind)
    start = t.segment_index(ind(n));
    if n==numel(ind)
        try
            stop =  t.segment_index(ind(n)+1);
        catch
            stop = t.no_snapshots;
        end
    else
        stop =  t.segment_index(ind(n)+1);
    end
    los = [los start:stop];
end

power = reshape( 10*log10( squeeze(sum( abs(cn.coeff).^2 , 3 )) ) , 4,[]);

mi = min(reshape(power,[],1)) - 5;
ma = max(reshape(power,[],1)) + 5;

ar = ones(1,cn.no_snap) * ma;
ar(los) = mi;

figure('Position',[ 100 , 100 , 1000 , 700]);
a = area(dist,ar);
set(a(1),'FaceColor',[0.7 0.9 0.7]);
set(a,'LineStyle','none')

hold on
plot(dist,power')
hold off

xlabel('Track [m]');
ylabel('Received Power per MIMO LINK [dB]');
axis([0 t.get_length mi ma])
legend('LOS','P_{11}','P_{12}','P_{21}','P_{22}',4)
box on

title('Received power along the track')

%%
% The next plot shows the RMS delay spread along the path for the first
% MIMO element. Again, shaded ares are for the LOS segments.

pow_tap = abs(squeeze(cn.coeff(1,1,:,:))).^2;
pow_sum = sum( pow_tap,1 );
mean_delay = sum( pow_tap.*cn.delay ,1) ./ pow_sum;
ds = sqrt( sum( pow_tap.*cn.delay.^2 ,1)./ pow_sum - mean_delay.^2 );
ar = zeros(1,cn.no_snap);
ar(los) = 10000;

figure('Position',[ 100 , 100 , 1000 , 700]);
a = area(dist,ar);
set(a(1),'FaceColor',[0.7 0.9 0.7]);
set(a,'LineStyle','none')

hold on
plot( dist , ds*1e9  )
hold off

ma = 1e9*( max(ds)+0.1*max(ds) );

axis([0 t.get_length 0 ma])
xlabel('Track [m]');
ylabel('Delay Spread [ns]');
legend('LOS','\sigma_\tau',1)
title('Position dependant delay spread');



%%
disp(['QuaDRiGa Version: ',simulation_parameters.version])
