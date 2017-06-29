function generate_initial_paths(h_channel_builder)
%GENERATE_INITIAL_PATHS Generate delays and powers (private)
%
% The generation of initial delays is done according to the winner WIM2
% model. The only difference is that here we assume that the LOS component
% is always present. In a heavily shaded environment, however, the K-factor
% becomes very small (below 0). Delays are drawn randomly from a delay
% distribution defined in the parameter set. This distribution can be
% estimated from measurements and is scenario dependent.
%
% Since a delay represents a longer path-distance, the signal will also
% have a higher attenuation. The NLOS cluster powers are calculated
% assuming a single slope exponential power delay profile. Power assignment
% depends on the delay spread and the per cluster shadow fading which is
% predefined in the parameter database.
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

use_ground_reflection = h_channel_builder.par.simpar.use_ground_reflection;
if use_ground_reflection
    L = h_channel_builder.NumClusters - 1;                  % no. taps
else
    L = h_channel_builder.NumClusters;
end
N = h_channel_builder.par.no_positions;                     % no. positions
ds = h_channel_builder.par.ds.';

if L ~= 1                                                   % LOS and NLOS components
    
    % First, we create the path delays
    taus = ds(:, ones(1, L));                               % initialize output with RMS-DS
    taus = -h_channel_builder.par.scenpar.r_DS * taus .* log(rand(N, L)); % generate random log-normal variables
    taus = sort(taus, 2);                                  	% sort output
    taus = taus - taus(:, ones(1, L));                      % normalize min. delay to zero
    
    % Generate the cluster powers
    % Note that the K-Factor correction is not yet applied to the delays.
    
    pow = exp(-taus .*...
        (h_channel_builder.par.scenpar.r_DS - 1)  ./...
        (h_channel_builder.par.scenpar.r_DS * ds(:, ones(1, L)))  );
    
    % Apply per-path shadowing
    pow = pow .* 10.^(-randn(N, L) * h_channel_builder.par.scenpar.LNS_ksi * 0.1);
    h_channel_builder.pow_wo_kf = pow;
    
    % Apply K-Factor in formation to cluster powers.
    pow(:, 1) = h_channel_builder.par.kf' .* sum(pow(:, 2:L), 2);
    
    % Normalize the path power such that the sum is equal to 1
    for n = 1 : N
        pow(n, :) = pow(n, :) / sum(pow(n, :));
    end
    
    % Normalize the taus to match the given delay spread
    ds_out = abs(sqrt(sum(pow.*(taus.^2), 2)  - sum(pow.*taus, 2).^2)  );
    scale = ds./ds_out;
    taus = scale(:, ones(1, L)).*taus;
    
else
    pow = ones(N, 1);
    taus = zeros(N, 1);
end

if h_channel_builder.par.simpar.use_ground_reflection
    
    % Calculate the delay of the ground reflection
    txpos = h_channel_builder.par.tx_position;
    rxpos = h_channel_builder.par.positions;
    gfpos = [ rxpos([1,2],:) ; -rxpos(3,:) ];
    
    % EoD for the GR
    theta_gr = -atan( ( -rxpos(3,:) - txpos(3) ) ./ ...
        sqrt( (txpos(1) - rxpos(1,:)).^2 +...
        (txpos(2) - rxpos(2,:)).^2 ) );
    
    oN = ones(1,N);
    d_3d = sqrt( sum((rxpos - txpos(:,oN)).^2,1) ).';
    d_gf = sqrt( sum((gfpos - txpos(:,oN)).^2,1) ).';
    tau_gr = ( d_gf-d_3d ) / h_channel_builder.par.simpar.speed_of_light;
       
    % Calculate the scaling facor for the ground reflection
    ii = ( taus - tau_gr(:,ones(L,1)) ) < 0;
    
    % The sum-power of the NLOS paths arriving before the ground reflection
    tmp = pow;
    tmp(~ii) = 0;
    tmp(:,1) = 0;
    P_before_gr  = sum(tmp,2);
    
    % The reflection coefficient
    epsilon_r = 15;                             % Dielectric Constant
    Zh = sqrt(epsilon_r - cos(theta_gr).^2);    % Horizintal Polarization
    Zv = Zh./epsilon_r;                         % Vertical Polarization
    R = 0.5*( sin( theta_gr ) - Zv ) ./ ( sin( theta_gr ) + Zv ) +...
        0.5*( sin( theta_gr ) - Zh ) ./ ( sin( theta_gr ) + Zh );
    R = abs(R);
    
    x = R/(1+R)*pow(:,1);
    %x = (R/(1+R)*pow(:,1) - P_before_gr)./pow(:,1);
    x(x<0) = 0;
    
    % Calculate the delay scaling factor
    a = ds.^2;
    b = pow(:,1) .* tau_gr.^2;
    c = sum( pow(:,2:end) .* taus(:,2:end).^2 , 2 );
    d = pow(:,1) .* tau_gr;
    e = sum( pow(:,2:end) .* taus(:,2:end) , 2 );
    y = ( sqrt( a.*c - a.*e.^2 - b.*c.*x + b.*e.^2.*x + c.*d.^2.*x.^2 ) + d.*e.*x ) ./ (c-e.^2);
    
    ii = imag(y)~=0;
    x( ii ) = 0;
    y( ii ) = 1;
    
    % Update the powers and dleays including the ground reflection
    pow  = [ pow(:,1).*(1-x) , pow(:,1).*x  , pow(:,2:end)                     ];
    taus = [ taus(:,1)       , tau_gr       , taus(:,2:end).*y(:,ones(1,L-1))  ];
end  

h_channel_builder.taus = taus;
h_channel_builder.pow = pow;

end
