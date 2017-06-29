%% Drifting Phases and Delays
%
% Drifting is the method used for obtaining time evolution within one
% state. This tutorial demonstrates the effect of â€œdriftingâ€?on the channel
% coefficients. It shows how drifting can be enabled and disabled as well
% as how the resulting data can be analyzed.
%
% Drifting is an essential feature of the channel model. Drifting
% enables a continuous time evolution of the path delays, the path phases,
% the departure- and arrival angles and the LSPs. It is thus the enabling
% feature for time continuous channel simulations. Although drifting was
% already available in the SCME branch of the WINNER channel model, it did
% not make it into the main branch. Thus, drifting is not available in the
% WIM1, WIM2 or WIM+ model. Here the functionality is implemented again.
% This script focuses on the delay and the phase component of the drifting
% functionality.        


%% Channel model setup and coefficient generation
% First, we parameterize the channel model. We start with the basic
% simulation parameters. For the desired output, we need two additional
% options: we want to evaluate absolute delays and we need to get all 20
% subpaths. Normally, the subpaths are added already in the channel
% builder.

clear all
close all

set(0,'defaultTextFontSize', 14)
set(0,'defaultAxesFontSize', 14)

s = simulation_parameters;
s.center_frequency = 2.53e9;
s.sample_density = 4;
s.map_resolution = 2;
s.use_absolute_delays = 1;

%%
% Second, we define a user track. Here we choose a linear track with a
% length of 30 m. The track start 20 m east of the transmitter and runs in
% east direction, thus linearly increasing the distance from the receiver.

l = layout( s );
l.tx_position(3) = 25;
l.track.generate('linear',30,0);
l.track.initial_position = [20;0;0];
l.track.scenario = 'WINNER_UMa_C2_LOS';
l.track.interpolate_positions( s.samples_per_meter );
l.visualize;
view(-33,45);

%%
% Now, we generate the LSPs. To get repeatable results, we set a specific
% random seed. This is a MATLAB internal function and is not a feature of
% the channel model. We also set the shadow fading and K-factor to 1 and
% disable the path loss model.

RandStream.setGlobalStream(RandStream('mt19937ar','seed',5));
p = l.create_parameter_sets;
p.sf_map = 0;
p.kf_map = 0;
p.plpar = [];                               % Disable path loss model
p.update_parameters;


%%
% Now, we generate the channel coefficients. The first run uses the new
% drifting module, the second run disables it. Note that drifting needs
% significantly more computing resources. In some scenarios it might thus
% be useful to disable the feature to get quicker simulation results.

s.drifting_precision = 5;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',2));
c = p.get_channels;

s.drifting_precision = 0;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',2));
d = p.get_channels;


c.individual_delays = 0;

%% Results and discussion
% The following plots represent the results of the test.

figure
distance = 20+(1:c.no_snap)*l.track.get_length/c.no_snap;
plot( distance, c.delay(1,:)*1e9 , '-b' )
hold on
plot( distance, d.delay(1,:)*1e9 , '-.b' )
plot( distance, c.delay(2,:)*1e9 , '-r' )
plot( distance, d.delay(2,:)*1e9 , '-.r' )
hold off
xlabel('Distance from track start point')
ylabel('Delay [ns] ')


%%
% The first plot shows the delay of the LOS tap (blue) and the delay of the
% first NLOS tap (red) vs. distance. The solid lines are from the channel
% with drifting, the dashed lines are from the channel without. The LOS
% delay is always increasing since the Rx is moving away from the Tx.
% However, the increase is not linear due to the 25 m height of the Tx.
% Without drifting, the delays are not updated and stay constant during the
% segment. The position of the first scatterer is in close distance to the
% Rx (only some m away). When moving, the Rx first approaches the scatterer
% (delay gets a bit smaller) and then the distance increases again.
figure
phase = unwrap(angle(squeeze((c.coeff(1,1,2,:,:)))));
plot( distance,phase )
xlabel('Distance from track start point')
ylabel('Continuous phase')

%%
% The second plot shows the phases of the 20 subpaths of the first NLOS tap
% for the drifting case. Note that the phases are not linear. This comes
% from the close proximity of the scatterer to the initial Rx position. The
% position of all 20 reflection points are calculated by the channel model.
% Those position mark the position of the last bounce scatterer (LBS). When
% moving the Rx, the distance to the LBS changes for each subpath and so
% does the phase. Here, the phase of each of the subpaths is calculated
% from the length of the path.
figure
pow = abs(squeeze(sum( c.coeff(1,1,2,:,:) , 5 ))).^2;
plot( distance,10*log10(pow),'-r' )
xlabel('Distance from track start point')
ylabel('Tap power (dB)')


%%
% This plot shows the power of the first NLOS tap along the track. The
% fading is significantly higher in the beginning and becomes much less
% strong towards the end.
figure
phase = unwrap(angle(squeeze((d.coeff(1,1,2,:,:)))));
plot( distance,phase )
xlabel('Distance from track start point')
ylabel('Continuous phase')

%%
% Without drifting, the phases of the subpaths are approximated by assuming
% that the angles to the LBSs do not change. However, this only holds when
% the distance to the LBS is large. Here, the initial distance is small
% (ca. 5 m). When the initial angles are kept fixed along the track, the
% error is significant. Here, the phase ramp is negative, indicating a
% movement direction towards the scatterer and thus a higher Doppler
% frequency. However, when the scatterer is passed, the Rx moves away from
% the scatterer and the Doppler frequency becomes lower. This is not
% reflected when drifting is turned off.
%
% Note here, that with shorter delay spreads (as e.g. in satellite
% channels), the scatterers are placed closer to the Rxs initial position.
% This will amplify this effect. Hence, for correct time evolution results,
% drifting needs to be turned on.


%%
figure
pow = abs(squeeze(sum( d.coeff(1,1,2,:,:) , 5 ))).^2;
plot( distance,10*log10(pow),'-r' )
xlabel('Distance from track start point')
ylabel('Tap power (dB)')

%%
%close all
disp(['QuaDRiGa Version: ',simulation_parameters.version])
