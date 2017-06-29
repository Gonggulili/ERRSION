%% Visualizing RHCP/LHCP Patterns
%
% The internal algorithms of the channel model only work with linear
% polarization. The antenna patterns are thus only stored in H/V
% polarization. This tutorial defines two circular patch antennas and places
% them in an environment. It then rotates the transmit antenna degree by
% degree and thus samples all azimuth and elevation angles. The channel
% model is set up to record the channel response and thus record the
% RHCP/LHCP response like in a measurement in an anechoic chamber.

%% Set up the array
% We set up a patch antenna with an opening angle of 90°. We then copy
% that patch and rotate it by 90° around the x-axis to create an X-Pol
% array and point the entire array skywards. We then set the coupling to
% +/- 90° phase to transmit circular polarized waves.

close all
clear all

resolution = 10;                        % in Degrees

a = array('custom',90,90,0);
a.set_grid( (-180:resolution:180)*pi/180 , (-90:resolution:90)*pi/180 );
a.copy_element(1,2);
a.rotate_pattern(90,'x',2);
a.coupling = 1/sqrt(2) * [1 1;1j -1j];  % LHCP and RHCP

b = a.copy_objects;                     % Circular receiver
b.field_pattern_vertical(:,:,1) = 1;
b.field_pattern_vertical(:,:,2) = 0;
b.field_pattern_horizontal(:,:,1) = 0;
b.field_pattern_horizontal(:,:,2) = 1;

a.rotate_pattern(-90,'y');              % Point transmitter skywards


%% Place arrays in layout
% We place two of those arrays in a layout. The scenario 'LOSonly' has no
% NLOS scattering. One can see this setup as a perfect anechoic chamber.

l = layout;
l.simpar.show_progress_bars = 0;
l.simpar.drifting_precision = 0;

l.rx_position = [10;0;0];
l.tx_position = [0;0;0];
l.track.no_snapshots = 1;
l.track.ground_direction = pi;
l.track.scenario = 'LOSonly';
l.tx_array = a;
l.rx_array = b;

p = l.create_parameter_sets;
[~,cb] = p.get_channels;
cb.pin = zeros( size(cb.pin) );

%% Get array response
% We now sample the array response for each degree in the antenna array.

pat = zeros( a.no_el , a.no_az , 2 , 2 );

values = a.no_az;
fprintf('Calculating  ['); m0=0; tStart = clock;  % A Status message
for n = 1:a.no_az
    m1=ceil(n/values*50); if m1>m0; for m2=1:m1-m0; fprintf('o'); end; m0=m1; end;
    
    a1 = a.copy_objects;
    a1.rotate_pattern( a.azimuth_grid(n)*180/pi , 'z');
    
    for m=1:a.no_el
        a2 = a1.copy_objects;
        a2.rotate_pattern( a.elevation_grid(m)*180/pi,'y');
        cb.par.tx_array = a2;
        c = cb.get_channels;
        pat(m,n,:,:) = c.coeff;
    end
end
fprintf('] %5.0f seconds\n',round( etime(clock, tStart) ));

%% Plot
% For plotting we use the internal function of the array class. We adjust
% the title of the figures accordingly.

d = a.copy_objects;
d.field_pattern_vertical =  pat(:,:,:,1) ;
d.field_pattern_horizontal = pat(:,:,:,2) ;
d.visualize


%%
disp(['QuaDRiGa Version: ',simulation_parameters.version])
