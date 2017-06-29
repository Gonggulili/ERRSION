%% Geometric Polarization
%
% This tutorial shows how to study polarization effects with QuaDRiGa.
% Different linearly polarized antennas are defined at the transmitter and
% the receiver, the channel between them is calculated and the polarization
% effects are evaluated.
%
% We demonstrate the polarization rotation model that calculates the path
% power for polarized antenna arrays. We do this by setting up the
% simulation with different H/V polarized antennas at the transmitter and
% at the receiver. Then we define a circular track around the receiver.
% When the receiver moves around the transmitter, it changes its antenna
% orientation according to the movement direction. In this way, all
% possible departure and elevation angles are sampled. Depending on the
% antenna orientation, the polarizations are either aligned (e.g. the Tx is
% V-polarized and the Rx is V-polarized), they are crossed (e.g. the Tx is
% V-polarized and the Rx is H-polarized), or the polarization orientation
% is in between those two. The generated channel coefficients should
% reflect this behavior.

%% Setting up the simulation environment
% First, we have to set up the simulator with some default settings. Here,
% we choose a center frequency of 2.1 GHz. We also want to use drifting in
% order to get the correct angles for the LOS component and we set the
% number of transmitters and receivers to one.

close all
clear all

set(0,'defaultTextFontSize', 14)
set(0,'defaultAxesFontSize', 14)

s = simulation_parameters;                      % Set the simulation parameters
s.center_frequency = 2.1e9;                     % Center-frequency: 2.1 GHz
s.use_polarization_rotation = 1;
s.samples_per_meter = 360/(40*pi);              % One sample per degree
s.drifting_precision = 1;


%% Setting up the antenna arrays
% In the second step, we set up our antenna arrays. We use the synthetic
% dipole antennas for this case. Those antennas show perfect polarization
% characteristics. First, we generate a single dipole with V-polarization.
% Then, we create multiple copies of this antenna element and rotate them
% by 45 and 90 degrees, respectively. We then use the same array for the
% receiver.

l = layout( s );                                % Create a new Layout
l.tx_array.generate('dipole');                  % create V-polarized dipole
l.tx_array.set_grid( (-180:10:180)*pi/180 , (-90:10:90)*pi/180 );

l.tx_array.field_pattern_vertical = ...         % Normalize
    l.tx_array.field_pattern_vertical / max(max(l.tx_array.field_pattern_vertical));

l.tx_array.set_grid( (-180:5:180)*pi/180 , (-90:5:90)*pi/180 );
l.tx_array.copy_element(1,2:3);                 % Duplicate the element two more times
l.tx_array.rotate_pattern(45,'y',2);            % 45deg polarization
l.tx_array.rotate_pattern(90,'y',3);            % 90deg polarization
l.tx_array.visualize;                           % Plot the array
l.rx_array = l.tx_array;                        % Use the same array for the Rx


%% Defining a track
% The third step defines the track. Here, we use a circle with 40 m diameter
% starting in the east, traveling north. We also choose a LOS scenario
% since we want to study the LOS polarization. The transmitter is located
% 12 m north of the center of the circle at an elevation of 6 m.

l.tx_position = [ 0 ; 12 ; 6 ];                 % Tx position
l.rx_position = [ 20 ; 0 ; 0 ];                 % Start position for the Rx track
l.track.generate('circular',40*pi,0);           % A circular track with radius 20 m
l.track.scenario = 'BERLIN_UMa_LOS';            % Chosse the Urban Macro LOS scenario
l.track.interpolate_positions( s.samples_per_meter );     % Interpolate positions
l.visualize;                                    % Plot the track


%% Generating channel coefficients
% Now, we have finished the parametrization of the simulation and we can
% generate the parameters. We thus create a new set of correlated LSPs and
% the fix the shadow fading and the K-factor to some default values. This
% disables the drifting for those parameters. We need to do that since
% otherwise, drifting and polarization would interfere with each other.

RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));

p = l.create_parameter_sets;                    % Create parameter sets
p.kf_map = 3;                                   % Fix KF to 3 dB
p.sf_map = 0;                                   % Fix SF to 0 dB
p.plpar = [];                                   % Disable path loss model
p.update_parameters;

[c,cb] = p.get_channels;                        % Get the channel coefficients

%% Results and Evaluation
% We now check the results and confirm, if they are plausible or not. We
% start with the two vertically polarized dipoles at the Tx and at the Rx
% side.
%
% The model creates 15 taps, which is the default for the BERLIN_UMa_LOS
% scenario. Without path-loss and shadow fading (SF=1), the power is
% normalized such that the sum over all taps is 1W and with a K-Factor of 3
% dB, we get a received power of 0.67W for the LOS component. The remaining
% 0.33W are in the NLOS components. The results can be seen in the
% following figure.

figure                                          % New figure
plot(abs(squeeze( c.coeff(1,1,:,:) )').^2);     % Plot the graph
axis([0 360 -0.1 1]);                           % Set the axis
xlabel('Position [degrees]');                   % Add description
ylabel('LOS Power, linear scale');
title('Tx: Vertical , Rx: Vertical');           % Add title

disp(['LOS power:  ',num2str(mean( abs(c.coeff(1,1,1,:)).^2 , 4))])
disp(['NLOS power: ',num2str(mean( sum(abs(c.coeff(1,1,2:end,:)).^2,3) , 4))])

%%
% The LOS power is almost constant when the Rx is south of the Tx. However,
% in close proximity (at 90deg), the power is lowered significantly. This
% comes from the 6 m elevation of the Tx. When the Rx is almost under the
% Tx, the radiated power of the Dipole is much smaller compared to the
% broadside direction. The average power of the LOS is thus also lowered to
% 0.56 W. The average sum-power if the 7 NLOS components is 0.26 W. This
% mainly come from the XPR which leakes some power from the vertical- into
% the horizontal polarization and thus reduces the received power on the
% vertically polarized Dipole.
%
% Next, we study two cases. Either the Tx is vertical polarized and the Rx
% is at 45deg or vise versa.

figure                                          % New figure
plot(abs(squeeze( c.coeff(2,1,1,:) )).^2);      % Tx vertical, Rx 45deg
hold on
plot(abs(squeeze( c.coeff(1,2,1,:) )).^2,'--r');	% Tx 45deg, Rx vertical
hold off
axis([0 360 -0.1 1]);
legend('Tx vertical, Rx 45\circ', 'Tx 45\circ, Rx vertical')
xlabel('Position [degrees]');
ylabel('LOS Power, linear scale');
title('Tx: Vertical , Rx: 45\circ');

%%
% The receiver changes its direction in a way that it always has the same
% orientation towards the Tx. However, due to the displacement of the Tx,
% the radiated power towards the Tx becomes minimal at around 90deg. This
% minimum is visible in both curves (blue and red). However, the pole of
% the 45deg slanted dipole now points to a different direction which explains
% the difference in the two lines. When the Rx is at 45deg and the Tx is
% vertical, the pole is in the right half if the circle - resulting in a
% lower received power. When the Rx is Vertical and the Tx is 45deg, the
% minimum power is achieved in the left half of the circle.
%
% Next, we evaluate the two dipoles which are rotated by 45deg. When moving
% around the circle, the Tx stays fixed and the Rx rotates. Subsequently,
% at one position, we will have both dipoles aligned and at another
% position, both will be crossed. When they are crossed, the received power
% will be 0 and when they are aligned, the power will match the first plot
% (two vertical dipoles). This can be seen in the following figure.

figure                                          % New figure
plot(abs(squeeze( c.coeff(2,2,1,:) )).^2 , 'Linewidth',1);
axis([0 360 -0.1 1]);
set(gca,'XTick',0:45:360)
xlabel('Position on circle [degrees]');
ylabel('LOS Power, linear scale');
title('Tx: 45\circ , Rx: 45\circ');


%%
% In the last figure, we have the Tx-antenna turned by 90deg. It is thus
% lying on the side and it is horizontally polarized. For the Rx, we
% consider three setups: Vertical (blue line), 45deg (green line) and 90deg
% (red line). Note that the Tx is rotated around the y-axis. At the initial
% position (0deg), the Rx (45 and 90deg) is rotated around the x-axis. This is
% because the movement direction.

figure                                          % New figure
plot(abs(squeeze( c.coeff(:,3,1,:) ))'.^2);
axis([0 360 -0.1 1]);
legend('Rx: 0\circ','Rx: 45\circ','Rx: 90\circ' )
xlabel('Position [degrees]');
ylabel('LOS Power, linear scale');
title('Tx: 90\circ , Rx: 0\circ, 45\circ, 90\circ');


%%
% When the receiver is vertical (blue line), both antennas are always
% crossed. There is no position around the circle where a good link can be
% established. When the receiver is horizontal (red line), however, there
% are two points where the two dipoles are aligned. For the 45deg dipole, the
% same behavior can be observed but with roughly half the power.

%%
%close all
disp(['QuaDRiGa Version: ',simulation_parameters.version])