%% Time Evolution and Scenario Transitions
%
% This tutorial shows how user trajectories, segments, and scenarios are
% defined. Channel coefficients are created for each segment separately.
% The channel merger combines these output into a longer sequence. The
% output sequences are evaluated for different settings of the model.
%
% The channel model generates the coefficients separately for each segment.
% In order to get a time-continuous output, these coefficients have to be
% combined. This is a feature which is originally described in the
% documentation of the WIM2 channel model, but which was never implemented.
% Since this component is needed for time-continuous simulations, it was
% implemented here. This script sets up the simulation and creates such
% time-continuous CIRs.

%% Channel model setup and coefficient generation
% First, we set up the channel model.
close all
clear all

set(0,'defaultTextFontSize', 14)
set(0,'defaultAxesFontSize', 14)

s = simulation_parameters;
s.center_frequency = 2.53e9;
s.sample_density = 4;
s.use_absolute_delays = 1;


%%
% Second, we create a more complex network layout featuring an elevated
% transmitter (25 m) and two receivers at 1.5 m height. The first Rx moves
% along a circular track around the receiver. The second receiver moves
% away from the Tx. Both start at the same point.
%
% Note here, that each track is split into three segments. The first Rx
% goes from an LOS area to a shaded area and back. The second
% track also start in the LOS area. Here, the scenario changes to another
% LOS segment and then to an NLOS segment. The LOS-LOS change will create
% new small-scale fading parameters, but the large scale parameters (LSPs)
% will be highly correlated between those two segments.

l = layout( s );                                % New layout
l.no_rx = 2;                                    % Two receivers

l.tx_array.generate('dipole');                  % Dipole antennas at all Rx and Tx
l.rx_array = l.tx_array;

l.tx_position(3) = 25;                          % Elevate Tx to 25 m

UMal = 'BERLIN_UMa_LOS';
UMan = 'BERLIN_UMa_NLOS';

l.track(1).generate('circular',20*pi,0);        % Circular track with 10m radius
l.track(1).initial_position    = [10;0;1.5];    % Start east, running nord
l.track(1).segment_index       = [1,40,90];     % Segments
l.track(1).scenario            = { UMal , UMan , UMal };

l.track(2).generate('linear',20,pi/8);          % Linear track, 20 m length
l.track(2).initial_position    = [10;0;1.5];    % Same start point
l.track(2).interpolate_positions( 128/20 );
l.track(2).segment_index       = [1,40,90];
l.track(2).scenario            = { UMal , UMal , UMan };

l.visualize;                                    % Plot all tracks
view(-33,45)

l.track.interpolate_positions( s.samples_per_meter );
l.track.compute_directions;


%%
% Now we create the channel coefficients. The fixing the random seed
% guarantees repeatable results (i.e. the taps will be at the same
% positions for both runs). Also note the significantly longer computing
% time when drifting is enabled.

RandStream.setGlobalStream(RandStream('mt19937ar','seed',2));
p = l.create_parameter_sets;

disp('Drifting enabled:');
s.drifting_precision = 1;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
c = p.get_channels;
cn = c.merge;

disp('Drifting disabled:');
s.drifting_precision = 0;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
d = p.get_channels;
dn = d.merge;


%% Results and discussion
% Now we plot the and discuss the results. We start with the power of the
% LOS tap along the circular track and compare the outcome with and without
% drifting.

degrees = (0:cn(1).no_snap-1)/cn(1).no_snap * 360;
los_pwr_drift = 10*log10(squeeze(abs(cn(1).coeff(1,1,1,:))).^2);
los_pwr_nodrift = 10*log10(squeeze(abs(dn(1).coeff(1,1,1,:))).^2);

figure
plot( degrees,los_pwr_drift )
hold on
plot(degrees,los_pwr_nodrift ,'-.r')
hold off

a = axis;
axis( [0 360 a(3:4) ] );

xlabel('Position on circle [degree]');
ylabel('Power of the LOS component');
title('Power of the LOS component for the circular track');
legend('time-continuous','no time-continuous','location','SouthEast');


%%
% When drifting is enabled (blue curve), the channel output after merging
% is time-continuous. The variations along the track come from the drifting
% K-Factor and the drifting shadow fading. When drifting is disabled, these
% parameters are not updated and kept fixed at their initial value.
%
% At the and of each segment, both channels are cross-faded. I.e. the power
% of the output of the first segment ramps down and the power of the second
% segment ramps up. Since drifting guarantees a time-continuous evolution
% of the phase, this ramping process is also time continuous and no
% artifacts are visible in the blue curve.
%
% Without drifting, however, the phases are approximated based on their
% initial values, the initial arrival- and departure angles and the
% traveled distance from the start point. However, since the Rx moves along
% a circular track, the angles change continuously which is not correctly
% modeled. The phase at the end of the first segment does not match the
% phase at the beginning of the second. When adding both components,
% artifacts appear as can be seen in the red curve.
%
% Next, we plot the power-delay profiles for both tracks. We calculate the
% frequency response of the channel and transform it back to time domain by
% an IFFT. Then, we create a 2D image of the received power at each
% position of the track. We start with the circular track.   

h = cn(1).fr( 100e6,512 );
h = squeeze(h);
pdp = 10*log10(abs(ifft(h,[],1).').^2);

figure
imagesc(pdp(:,1:256));
caxis([ max(max(pdp))-50 max(max(pdp))-5 ]);
colorbar

cm = colormap('hot');
colormap(cm(end:-1:1,:));

set(gca,'XTick',1:32:255);
set(gca,'XTickLabel',(0:32:256)/100e6*1e6);
xlabel('Delay [\mus]');

set(gca,'YTick',1:cn(1).no_snap/8:cn(1).no_snap);
set(gca,'YTickLabel', (0:cn(1).no_snap/8:cn(1).no_snap)/cn(1).no_snap * 360 );
ylabel('Position on circle [degree]');

title('PDP for the circular track with drifting');


%%
% The X-axis shows the delay in microseconds and the Y-axis shows the
% position on the circle. For easier navigation, the position is given in
% degrees. 0° means east (starting point), 90° means north, 180° west and
% 270° south. The LOS delay stays constant since the distance to the Tx is
% also constant. However, the power of the LOS changes according to the
% scenario. Also note, that the NLOS segment has significantly more taps
% due to the longer delay spread.      
%
% Next, we create the same plot for the linear track. Note the slight
% increase in the LOS delay and the high similarity of the first two LOS
% segments due to the correlated LSPs. Segment change is at around 6 m.

h = cn(2).fr( 100e6,512 );
h = squeeze(h);
pdp = 10*log10(abs(ifft(h,[],1).').^2);

figure
imagesc(pdp(:,1:256));
caxis([ max(max(pdp))-50 max(max(pdp))-5 ])
colorbar

cm = colormap('hot');
colormap(cm(end:-1:1,:));

set(gca,'XTick',1:32:255);
set(gca,'XTickLabel',(0:32:256)/100e6*1e6);
xlabel('Delay [\mus]');

set(gca,'YTick',1:cn(2).no_snap/8:cn(2).no_snap);
set(gca,'YTickLabel', (0:cn(2).no_snap/8:cn(2).no_snap)/cn(2).no_snap * 20 );
ylabel('Distance from start point [m]');

title('PDP for the linear track with drifting');


%%
% Last, we plot the same results for the linear track without drifting.
% Note here, that the LOS delay is not smooth during segment change. There
% are two jumps at 6 m and again at 13.5 m.

h = dn(2).fr( 100e6,512 );
h = squeeze(h);
pdp = 10*log10(abs(ifft(h,[],1).').^2);

figure
imagesc(pdp(:,1:256));
caxis([ max(max(pdp))-50 max(max(pdp))-5 ])
colorbar

cm = colormap('hot');
colormap(cm(end:-1:1,:));

set(gca,'XTick',1:32:255);
set(gca,'XTickLabel',(0:32:256)/100e6*1e6);
xlabel('Delay [\mus]');

set(gca,'YTick',1:cn(2).no_snap/8:cn(2).no_snap);
set(gca,'YTickLabel', (0:cn(2).no_snap/8:cn(2).no_snap)/cn(2).no_snap * 20 );
ylabel('Distance from start point [m]');

title('PDP for the linear track without drifting');


%%
%close all
disp(['QuaDRiGa Version: ',simulation_parameters.version])

%%
% The following code creates a movie of the time varying frequency response
% of the circular track. 



% h = cn(1).fr( 20e6,128 );
% h = squeeze(h);
% 
% mi = -90; ma = -80;
% while true
%     for n = 1:size(h,2)
%         pdp  = 10*log10(abs(h(:,n)).^2); 
%         plot(pdp)
%         ma = max( ma,max([pdp]) );
%         mi = min( mi,min([pdp]) );
%         axis([1,128,mi,ma])
%         title(n)
%         drawnow
%     end
% end

