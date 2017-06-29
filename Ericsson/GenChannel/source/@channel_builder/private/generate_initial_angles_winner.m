function generate_initial_angles_winner(h_channel_builder)
%GENERATE INITIAL_ANGLES Generate angular parameters (private)
%
% In three dimensional space, the signals originating from a cluster would
% arrive at the receiver under a certain angle in azimuth and elevation
% direction. The same hold for the angles of departure.
%
% Note: All generated angles are in radians!
%
% Reference: Kyösti, P.; Meinilä, J.; Hentilä, L. & others; IST-4-027756
% WINNER II D1.1.2 v.1.1: WINNER II Channel Models; 2007 
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if h_channel_builder.par.simpar.use_ground_reflection
    error('Ground reflection is not availabe for the WINNER angle mapping option');
end

N = h_channel_builder.par.no_positions;                   % no. positions
L = h_channel_builder.NumClusters;                        % no. taps
K = 10*log10(h_channel_builder.par.kf);              	  % Initial K-Factors in dB

% Get the correction values for the K-Factor and No-Clusters.
Cc = correction_function(L, K);

angles = h_channel_builder.par.get_angles*pi/180;         % Angles between Tx and Rx
pow = h_channel_builder.pow;                              % Path powers
asD = h_channel_builder.par.asD;                          % Azimuth Spread of Departure
esD = h_channel_builder.par.esD;                          % Elevation Spread of Departure
asA = h_channel_builder.par.asA;                          % Azimuth Spread of Arrival
esA = h_channel_builder.par.esA;                          % Elevation Spread of Arrival

AoD = zeros(N, L);                          % Azimuth departure angles
for n = 1:N
    pa = asD(n) / Cc(n) * 0.017453292519943;
    X = 2*randi(2, 1, L)-3;
    Y = randn(1, L) * asD(n)/7 * 0.017453292519943;
    
    AoD(n, :) = pa * sqrt(2) * sqrt(-log(pow(n, :)/max(pow(n, :))));
    AoD(n, :) = X.*AoD(n, :) + Y;
    AoD(n, 1) = 0; % LOS
    AoD(n, :) = AoD(n, :) + angles(1, n);
end
AoD = mod(AoD + pi, 2*pi) - pi;


EoD = zeros(N, L);                           % elevation departure angles
for n = 1:N
    pa = esD(n) / Cc(n) * 0.017453292519943;
    X = 2*randi(2, 1, L)-3;
    Y = randn(1, L) * esD(n)/7 * 0.017453292519943;
    
    EoD(n, :) = pa * sqrt(2) * sqrt(-log(pow(n, :)/max(pow(n, :))));
    EoD(n, :) = X.*EoD(n, :) + Y;
    EoD(n, 1) = 0; % LOS
    EoD(n, :) = EoD(n, :) + angles(3, n);
end
EoD = mod(EoD + pi, 2*pi) - pi;
EoD(EoD>pi/2) = pi - EoD(EoD>pi/2);
EoD(EoD<-pi/2) = -pi - EoD(EoD<-pi/2);


AoA = zeros(N, L);                           % azimuth arrival angles
for n = 1:N
    pa = asA(n) / Cc(n) * 0.017453292519943;
    X = 2*randi(2, 1, L)-3;
    Y = randn(1, L) * asA(n)/7 * 0.017453292519943;
    
    AoA(n, :) = pa * sqrt( -2*log(pow(n,:)/max(pow(n,:))) );
    AoA(n, :) = X.*AoA(n, :) + Y;
    AoA(n, 1) = 0; % LOS
    AoA(n, :) = AoA(n, :) + angles(2, n);
end
AoA = mod(AoA + pi, 2*pi) - pi;


EoA = zeros(N, L);                           % elevation arrival angles
for n = 1:N
    pa = esA(n) / Cc(n) * 0.017453292519943;
    X = 2*randi(2, 1, L)-3;
    Y = randn(1, L) * esA(n)/7 * 0.017453292519943;
    
    EoA(n, :) = pa * sqrt(2) * sqrt(-log(pow(n, :)/max(pow(n, :))));
    EoA(n, :) = X.*EoA(n, :) + Y;
    EoA(n, 1) = 0; % LOS
    EoA(n, :) = EoA(n, :) + angles(4, n);
end
EoA = mod(EoA + pi, 2*pi) - pi;
EoA(EoA>pi/2) = pi - EoA(EoA>pi/2);
EoA(EoA<-pi/2) = -pi - EoA(EoA<-pi/2);

h_channel_builder.AoD = AoD;
h_channel_builder.AoA = AoA;
h_channel_builder.EoD = EoD;
h_channel_builder.EoA = EoA;

end
