%% How to manually set SSPs in QuaDRiGa (Satellite Scenario)
%
% Sometimes it is necessary to adjust the small-scale fading directly. This
% tutorial creates a simple satellite scenario. Then, a list of large-scale
% parameters (LSPs) and their corresponding small-scale parameters (SSPs)
% is given. It is also demonstrated how those parameters are adjusted and
% how they influence the resulting channel coefficients.

%% Setting general parameters
% We set up some basic parameters such as center frequency and sample density.

close all
clear all

set(0,'defaultTextFontSize', 14)    % Set default font size for the plots
set(0,'defaultAxesFontSize', 14)

s = simulation_parameters;          % Basic simulation parameters
s.center_frequency      = 2.185e9;  % Center Frequency
s.sample_density        = 4;        % 4 samples / half wave length


%% Defining Antenna Arrays
% We set up our antenna arrays for the transmitter at the satellite and the
% receiver. We use synthetic patch antennas for this case. Two elements
% are crossed by an angle of 90 degree. The signal is then split and fed
% with a 90 degree phase shift to both elements generating RHCP and 
% LHCP signals.

% Create a patch antenna with 120 degree opening
a = array('custom',120,120,0);

% Copy element 1 to element 2.
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

%% Setting the satellite orbital position
% In this step, we set the antennas and the position of the
% satellite into a simulation layout. A layout object contains all the
% geometric information that are necessary to run the simulation. First,
% we define the position of the satellite. Since the model uses Cartesian
% coordinates, we have to transform the position of the satellite first.
%
% Choose a random satellite position (Astra 2, seen from Berlin).
% The distance only needs to be big enough to ensure insignificant changes
% in the reception angle on the ground.

sat_el      = 28.4;             % Elevation angle
sat_az      = 161.6;            % Azimuth angle (South = 180 degree)
rx_latitude = 51;               % Latitude of the Rx

% Approximate the satellite distance for GEO orbit
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

% Create a new layout
l = layout( s );

% Set the satellite position in the layout
l.tx_position = [ sat_x ; sat_y ; sat_z ];

% Set the transmit antenna in the layout
l.tx_array = a;

l.randomize_rx_positions( 500,0,0,0 );          % Random Rx position


%% Setting up a user track and assigning LSPs
% Here, we create a track with one segment. We then assign LSPs and
% automatically calculate SSPs. Those parameters can then be edited.

t = track('linear',5);                          % Linear track, 20 m length
t.interpolate_positions( s.samples_per_meter);  % Interpolate to sample density
t.initial_position = l.rx_position;
t.scenario = 'MIMOSA_10-45_LOS';                % Set scenarios

t.generate_parameters;                          % Generate LSPs

l.track = t;                    % Set track
l.rx_array = b;                 % Set Rx-Antenna

[~,~,cb] = l.get_channels;      % Get SSPs

%%
% The LSPs are stored in "t.par". The SSPs are stored in "cb". 
% The parameter are as follows:
% 
%    t.par.ds
%       The delay spread in [s]. This value is used to assign the powers
%       "cb.pow" and the delays "cb.taus" of each SC. You can only change
%       the value of "t.par.ds" before calling "[~,~,cb] = l.get_channels;".    
%
%    cb.taus
%       The individual delays in [s] for each SC.
%
%    t.par.kf
%       A vector of values for temporal evolution of the Ricean K-factor in
%       [dB] for each snapshot (sample point). The first value
%       "t.par.kf(1)" is used to scale the power of the first SC. The
%       vector is applied right before applying the antenna patterns.
%       Hence, it overwrites any changes you make to the SC powers in
%       "cb.pow". It should always be edited before scaling the individual
%       SC powers, if needed. 
% 
%    cb.pow
%       The power of each SC normalized to a sum-power of 1. The initial
%       Ricean K-factor is already applied. This is done in a way that the
%       following variables are identical:  
% 
%       cb.pow(1) / sum( cb.pow(2:end) )    % The KF in the SC powers
%       10^( 0.1*t.par.kf(1) )              % The automatically generated value
%       cb.par.kf                           % The initial value for the segment
%
%       Any changes you make will be overwritten by "t.par.kf". Hence, if
%       you edit "cb.pow", you must make sure to also correct the values in
%       "t.par.kf" and "cb.par.kf". This is done as follows:  

% Edit power values in cb.pow as you wish. Here, we scale the first and
% second SC.
cb.pow(1:2) = [0.5, 0.2];

% Make sure that the sum-power is 1
cb.pow = cb.pow / sum(cb.pow);

% Calcualte new KF
kf = 10*log10( cb.pow(1) / sum( cb.pow(2:end) ) );

% Correct KF in "t.par.kf" and "cb.par.kf"
t.par.kf = t.par.kf - ( t.par.kf(1) - kf );
cb.par.kf = 10.^( 0.1*t.par.kf(1) );

%%
%    t.par.pg
%       A vector of values for total path gain in [dB] for each snapshot.
%       This is equal to the sum-power (e.g. in [dBm]) minus the transmit
%       power (also in [dBm]). This vector is applied after the KF-scaling
%       and right before applying the antenna patterns.
%
%    t.par.asA
%       Azimuth angle spread of arrival in [deg]. You can only change
%       the value of "t.par.asA" before calling "[~,~,cb] = l.get_channels;".
%
%    cb.AoA
%       Azimuth arrival angles for each SC in global spherical coordinates
%       given in [rad]. The first angle corresponds to the satellite
%       azimuth position. The following values are identical:
%
%       -sat_az + 90            % Given in compass-coordinates
%       cb.AoA(1)/pi*180        % Given in global spherical coordinates
%
%    t.par.esA
%       Elevation angle spread of arrival in [deg]. You can only change
%       the value of "t.par.esA" before calling "[~,~,cb] = l.get_channels;".
%
%    cb.EoA
%       Elevation arrival angles for each SC in global spherical coordinates
%       given in [rad]. The first angle corresponds to the satellite
%       elevation position. The following values are identical:
%
%       sat_el
%       cb.EoA(1)/pi*180
%       
%    t.par.asD
%       Azimuth angle spread of departure in [deg]. This value is
%       irrelevant for satellite setups where the departure angle is almost
%       the same for each SC. You can only change the value of "t.par.asD"
%       before calling "[~,~,cb] = l.get_channels;".      
%
%    cb.AoD
%       Azimuth departure angles for each SC in global spherical coordinates
%       given in [rad]. All angles should point to the terminal position as
%       seen from the satellite. The following values are identical:
%
%       -sat_az + 90            % Given in compass-coordinates
%       cb.AoD/pi*180-180       % Given spherical coordinates
%
%    t.par.esD
%       Elevation angle spread of departure in [deg]. This value is
%       irrelavant for satellite setups where the departure angle is almost
%       the same for each SC.You can only change the value of "t.par.esA"
%       before calling "[~,~,cb] = l.get_channels;". 
%
%    cb.EoA
%       Elevation departure angles for each SC in global spherical
%       coordinates given in [rad]. All angles correspond roughly to the
%       negative satellite elevation position. The following values are
%       identical: 
%
%       -sat_el
%       cb.EoD/pi*180
%
% After making the required changes, you can combine the SC and the antenna
% patterns to calculate the channel coefficient by:

c = cb.get_channels;

%%
% A narrow band time-series is calculated by:

h = squeeze( sum(c.coeff,3) );

%%
% Plot for the power on each MIMO sub-channel.

power = abs( reshape( h,4,[] ) ).^2;
power = 10*log10(power).';

figure
[~,dist] = t.get_length;
plot(dist,t.par.pg,'--k','Linewidth',2)
hold on
plot(dist,power)
hold off
title('Simulated Power')
xlabel('Distance from start point [m]')
ylabel('Received Power [dBm]')
axis([0,5,min( power(:) )-5,max( power(:) )+2])
grid on
legend('LSP: t.par.pg','LHCP-LHCP','LHCP-RHCP','RHCP-LHCP','RHCP-RHCP',4);

%%
close all
disp(['QuaDRiGa Version: ',simulation_parameters.version])

