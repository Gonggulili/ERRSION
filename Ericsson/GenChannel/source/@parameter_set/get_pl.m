function [ loss , scale_sf ] = get_pl( h_parset , evaltrack , i_mobile )
%GET_PL Calculates the path loss
%
%   loss = GET_PL
%   calculates the path loss (PL) for all user positions. This function
%   implements the path loss models defined by the 'plpar' property of the
%   class 'parameter_set'.
%
%   [ loss , scale_sf ] = GET_PL
%   In some scenarios, the SF might change with increasing distance between
%   Tx and Rx. Hence, the shadow fading provided by the parameter map has
%   to be changed accordingly. The second output parameter "scale_sf" can
%   be used for scaling the (logarithmic) SF value from the map.
%
%   [ loss , scale_sf ] = GET_PL( evaltrack )
%   returns the PL along the positions of the provided track 'evaltrack'
%   instead of using the positions in the class property.
%
% !!!!   Modified version ebl:  !!!!!
%   Pathloss model "winner_los_h" added
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


if numel(h_parset) > 1
    error('??? Found array of "parameter_set".')
end

drifting_precision = h_parset.simpar.drifting_precision;
if nargin ~= 3
    drifting_precision = 0;
end

n_rx = 1;
n_tx = 1;
o_tx = 1;

txpos = h_parset.tx_position;

% Parse input variables
if nargin >= 2
    if ~( isa(evaltrack, 'track') )
        error('??? "evaltrack" must be of class "track".')
    end
    
    rxpos = [ evaltrack.positions(1,:) + evaltrack.initial_position(1) ; ...
        evaltrack.positions(2,:) + evaltrack.initial_position(2) ; ...
        evaltrack.positions(3,:) + evaltrack.initial_position(3) ];
    
    use_track = true;
    n_snapshots = evaltrack.no_snapshots;
    o_snapshots = ones(1,n_snapshots);
    
    if drifting_precision == 3
        e_rx = h_parset.rx_array(i_mobile).element_position;
        e_tx = h_parset.tx_array.element_position;
        
        n_rx = h_parset.rx_array(i_mobile).no_elements;
        n_tx = h_parset.tx_array.no_elements;
        
        o_rx = ones(1,n_rx);
        o_tx = ones(1,n_tx);
        
        gdir = evaltrack.ground_direction;
        c_gdir = cos(gdir);
        s_gdir = sin(gdir);
        
        e_rx_x = e_rx(1,:);
        e_rx_y = e_rx(2,:);
        e_rx_z = e_rx(3,:);
        
        % Apply the rotation
        q_s = zeros( 3,n_snapshots,n_rx );
        for i_rx = 1:n_rx
            q_s(1,:,i_rx) = c_gdir.*e_rx_x(i_rx) - s_gdir.*e_rx_y(i_rx);
            q_s(2,:,i_rx) = s_gdir.*e_rx_x(i_rx) + c_gdir.*e_rx_y(i_rx);
            q_s(3,:,i_rx) = e_rx_z(i_rx);
        end
        
        rxpos = q_s + rxpos(:,:,o_rx);
        txpos = reshape(e_tx,3,1,n_tx) + txpos(:,1,o_tx);
    end
else
    rxpos = h_parset.positions;
    use_track = false;
    n_snapshots = h_parset.no_positions;
    o_snapshots = ones(1,n_snapshots);
end

% Calculate the distance between Tx and Rx
d_3d = zeros(n_snapshots, n_rx, n_tx);
d_2d = zeros(n_snapshots, n_rx, n_tx);
for i_rx = 1:n_rx
    for i_tx = 1:n_tx
        d_3d(:, i_rx, i_tx) = ...
            sqrt(sum((rxpos(:, :, i_rx) - txpos(:, o_snapshots, i_tx)).^2 ));
        d_2d(:, i_rx, i_tx) = ...
            sqrt(sum((rxpos([1,2], :, i_rx) - txpos([1,2], o_snapshots, i_tx)).^2 ));
    end
end
d_3d = reshape(d_3d, 1, []);
d_2d = reshape(d_2d, 1, []);

% Calculate the Manhattan distance between Tx and Rx
d1_2d = zeros(n_snapshots, n_rx, n_tx);
d2_2d = zeros(n_snapshots, n_rx, n_tx);
for i_rx = 1:n_rx
    for i_tx = 1:n_tx
        d1_2d(:, i_rx, i_tx) = ...
            abs(rxpos(1, :, i_rx) - txpos(1, o_snapshots, i_tx)) ;
        d2_2d(:, i_rx, i_tx) = ...
            abs(rxpos(2, :, i_rx) - txpos(2, o_snapshots, i_tx));
    end
end
d1_2d = reshape(d1_2d, 1, []);
d2_2d = reshape(d2_2d, 1, []);

% Initialize output variables
loss = zeros(size(d_3d));
sf_sigma = loss;

% This implements the path loss models
if isfield( h_parset.plpar , 'model' )
    
    % Read commonly used data structures from class object (increases speed)
    par   = h_parset.plpar;
    CenterFrequency = h_parset.simpar.center_frequency/1e9;      % in GHz
    
    switch h_parset.plpar.model
        
        case 'logdist'
            loss        = par.A * log10(d_3d) + par.B + par.C * log10(CenterFrequency);
            
            if isfield( par , 'SF' )
                sf_sigma    = par.SF;
            else
                sf_sigma = h_parset.scenpar.SF_sigma + ...
                    h_parset.scenpar.SF_delta * ...
                    log10( h_parset.simpar.center_frequency / 1e9 );
            end
            
        case 'logdist_simple'
            loss        = par.A * log10(d_3d) + par.B;
            
            sf_sigma = h_parset.scenpar.SF_sigma + ...
                h_parset.scenpar.SF_delta * ...
                log10( h_parset.simpar.center_frequency / 1e9 );
            
        case 'constant'
            loss = par.A * ones( size(d_3d) );
            
            sf_sigma = h_parset.scenpar.SF_sigma + ...
                h_parset.scenpar.SF_delta * ...
                log10( h_parset.simpar.center_frequency / 1e9 );
            
        case 'winner_los'
            % From WINNER+ D5.3 pp. 74
            
            hBS = txpos(3);
            hMS = reshape( rxpos(3,:,:,o_tx) , 1 ,[] );
            
            hBS( hBS < 1.5  ) = 1.5;
            hMS( hMS < 1.5  ) = 1.5;
            
            % Calculate the breakpoint
            G = par.B1 + par.C1*log10(CenterFrequency) + par.D1*log10(hBS) + par.E1*log10(mean(hMS));
            H = par.B2 + par.C2*log10(CenterFrequency) + par.D2*log10(hBS) + par.E2*log10(mean(hMS));
            bp  = 10^( (H-G)/( par.A1-par.A2 ) );
            
            ind = d_2d<=bp;
            if sum(ind)>0
                loss(ind) = par.A1*log10(d_2d(ind)) + par.B1 + par.C1*log10(CenterFrequency)...
                    + par.D1*log10(hBS) + par.E1*log10(hMS(ind)) + par.F1*hMS(ind);
                sf_sigma(ind) = par.sig1;
            end
            
            ind = ~ind;
            if sum(ind)>0
                loss(ind) = par.A2*log10(d_2d(ind)) + par.B2 + par.C2*log10(CenterFrequency)...
                    + par.D2*log10(hBS) + par.E2*log10(hMS(ind)) + par.F2*hMS(ind);
                sf_sigma(ind) = par.sig2;
            end
            
        case 'winner_los_h'
            % From WINNER+ D5.3 pp. 74
            % Modification:  Full 3D distance is used !!!
            
            hBS = txpos(3);
            hMS = reshape( rxpos(3,:,:,o_tx) , 1 ,[] );
            
            hBS( hBS < 1.5  ) = 1.5;
            hMS( hMS < 1.5  ) = 1.5;
            
            % Calculate the breakpoint
            G = par.B1 + par.C1*log10(CenterFrequency) + par.D1*log10(hBS) + par.E1*log10(mean(hMS));
            H = par.B2 + par.C2*log10(CenterFrequency) + par.D2*log10(hBS) + par.E2*log10(mean(hMS));
            bp  = 10^( (H-G)/( par.A1-par.A2 ) );
            
            ind = d_3d<=bp;
            if sum(ind)>0
                pl(ind) = par.A1*log10(d_3d(ind)) + par.B1 + par.C1*log10(CenterFrequency)...
                    + par.D1*log10(hBS) + par.E1*log10(hMS(ind)) + par.F1*hMS(ind);
                sf_sigma(ind) = par.sig1;
            end
            
            ind = ~ind;
            if sum(ind)>0
                pl(ind) = par.A2*log10(d_3d(ind)) + par.B2 + par.C2*log10(CenterFrequency)...
                    + par.D2*log10(hBS) + par.E2*log10(hMS(ind)) + par.F2*hMS(ind);
                sf_sigma(ind) = par.sig2;
            end
            
        case 'winner_nlos'
            % From WINNER+ D5.3 pp. 74
            
            hBS = txpos(3);
            hMS = reshape( rxpos(3,:,:,o_tx) , 1 ,[] );
            
            hBS( hBS < 1.5  ) = 1.5;
            hMS( hMS < 1.5  ) = 1.5;
            
            if CenterFrequency < 1.5
                loss = ( par.A1 + par.Ah1 * log10( hBS ))*log10(d_2d) + par.B1 + ...
                    par.C1*log10(CenterFrequency) + ...
                    par.D1*log10(hBS) + ...
                    par.E1*log10(hMS) + ...
                    par.F1*hMS;
                
            elseif CenterFrequency >= 1.5 && CenterFrequency < 2
                loss = ( par.A2 + par.Ah2 * log10( hBS ))*log10(d_2d) + par.B2 + ...
                    par.C2*log10(CenterFrequency) + ...
                    par.D2*log10(hBS) + ...
                    par.E2*log10(hMS) + ...
                    par.F2*hMS;
                
            else % CenterFrequency >= 2
                loss = ( par.A3 + par.Ah3 * log10( hBS ))*log10(d_2d) + par.B3 + ...
                    par.C3*log10(CenterFrequency) + ...
                    par.D3*log10(hBS) + ...
                    par.E3*log10(hMS) + ...
                    par.F3*hMS;
            end
            sf_sigma = h_parset.scenpar.SF_sigma + ...
                h_parset.scenpar.SF_delta * ...
                log10( h_parset.simpar.center_frequency / 1e9 );
            
        case 'winner_pathloss'
            % See WINNER II D1.1.2 V1.2 (2007-09) p43 Equation (4.23)
            % PL/[dB] = A log10(d/[m]) + B + C log10(fc/[GHz]/5) + X
            
            loss = par.A * log10(d_3d) + par.B +...
                par.C * log10(h_parset.simpar.center_frequency/5e9) + par.X;
            
            if isfield( par , 'SF' )
                sf_sigma    = par.SF;
            else
                sf_sigma = h_parset.scenpar.SF_sigma + ...
                    h_parset.scenpar.SF_delta * ...
                    log10( h_parset.simpar.center_frequency / 1e9 );
            end
            
        case 'umts_vehicular'
            % See UMTS 30.03 version 3.1.0 (1997-11) p. 41
            % Section B.1.4.1.3 Path Loss Model for Vehicular Test
            % Environment
            % PL/[dB] = 40*(1-4*1e-3*(delta_hB/[m]))*log10(d/[m]/1000)-18*log10(delta_hB/[m]) + 21*log10(fc/[GHz]*1000) + 80;
            % assume height difference of 15 meters as defined in tr 36.814
            
            delta_hB = 15;
            
            loss = 40*(1-4*1e-3*(delta_hB))*log10(d_3d/1000)-18*log10(delta_hB) + 21*log10(CenterFrequency*1000) + 80;
            sf_sigma = h_parset.scenpar.SF_sigma + ...
                h_parset.scenpar.SF_delta * ...
                log10( h_parset.simpar.center_frequency / 1e9 );
            
        case 'itu_nlos'
            
            loss = 140.7+36.7*log10(d_3d*1e-3);
            sf_sigma = h_parset.scenpar.SF_sigma + ...
                h_parset.scenpar.SF_delta * ...
                log10( h_parset.simpar.center_frequency / 1e9 );
            
        case '3gpp_3d_uma_los'
            
            hBS = txpos(3);
            hMS = reshape(rxpos(3, :, :, o_tx), 1, []);
            
            hBS(hBS < 1.5) = 1.5;
            hMS(hMS < 1.5) = 1.5;
            
            % compute effective antenna heights. Assume that the
            % environment height hE equals to 1.0m.
            hE = 1;
            hBS_eff = hBS - hE;
            hMS_eff = hMS - hE;
            
            % calculate the break point (f_c in Hz!!!)
            c = h_parset.simpar.speed_of_light;
            bp = 4*(hBS_eff)*(hMS_eff)*CenterFrequency/c*1e9;
            
            ind = d_2d <= bp;
            if sum(ind) > 0 % for users < break point
                loss(ind) = par.A1*log10(d_3d(ind)) + ...
                    par.B1 + ...
                    par.C1*log10(CenterFrequency) + ...
                    par.D1*log10(bp(ind).^2 + (hBS - hMS(ind)).^2);
                
                sf_sigma(ind) = par.sig1;
            end
            
            ind = ~ind;
            if sum(ind) > 0 % for users > break point
                loss(ind) = par.A2*log10(d_3d(ind)) + ...
                    par.B2 + ...
                    par.C2*log10(CenterFrequency) + ...
                    par.D2*log10(bp(ind).^2 + (hBS - hMS(ind)).^2);
                
                sf_sigma(ind) = par.sig2;
            end
            
        case '3gpp_3d_uma_nlos'
            
            hBS = txpos(3);
            hMS = reshape(rxpos(3, :, :, o_tx), 1, []);
            
            hBS(hBS < 1.5) = 1.5;
            hMS(hMS < 1.5) = 1.5;
            
            % Calculate the breakpoint
            c = h_parset.simpar.speed_of_light;
            hE = 1; % m
            
            % break point; use f_c in Hz!!!
            bp = 4*(hBS-hE)*(hMS-hE)*CenterFrequency/c*1e9;
            
            ind = d_2d <= bp;
            
            loss_1 = zeros(size(ind));
            sf_sigma_tmp = loss_1;
            
            if sum(ind) > 0
                loss_1(ind) = par.A1*log10(d_3d(ind)) + ...
                    par.B1 + ...
                    par.C1*log10(CenterFrequency) + ...
                    par.D1*log10(bp(ind).^2 + (hBS - hMS(ind)).^2);
                sf_sigma_tmp(ind) = par.sig1;
            end
            
            ind = ~ind;
            if sum(ind) > 0
                loss_1(ind) = par.A2*log10(d_3d(ind)) + par.B2 + par.C2*log10(CenterFrequency)...
                    + par.D2*log10(bp(ind).^2 + (hBS - hMS(ind)).^2);
                sf_sigma_tmp(ind) = par.sig2;
            end
            
            loss_2 = par.B3 + ...
                - (24.37 - 3.7*(20/hBS).^2)*log10(hBS) ...
                + (par.A3 -3.1*log10(hBS))*(log10(d_3d) - 3) ...
                + par.C3*log10(CenterFrequency) ...
                - 0.6*(hMS - 1.5);
            
            use_loss_1 = bsxfun(@gt, loss_1, loss_2);
            loss(use_loss_1) = loss_1(use_loss_1);
            sf_sigma(use_loss_1) = sf_sigma_tmp(use_loss_1);
            
            use_loss_2 = ~use_loss_1;
            loss(use_loss_2) = loss_2(use_loss_2);
            sf_sigma(use_loss_2) = par.sig3;
            
        case '3gpp_3d_umi_nlos'
            
            hBS = txpos(3);
            hMS = reshape(rxpos(3, :, :, o_tx), 1, []);
            
            hBS(hBS < 1.5) = 1.5;
            hMS(hMS < 1.5) = 1.5;
            
            % Calculate the breakpoint
            c = h_parset.simpar.speed_of_light;
            
            hE = 1; % m
            % break point; use f_c in Hz!!!
            bp = 4*(hBS-hE)*(hMS-hE)*CenterFrequency/c*1e9;
            
            ind = d_2d <= bp;
            
            loss_1 = zeros(size(ind));
            sf_sigma_tmp = loss_1;
            
            if sum(ind) > 0
                loss_1(ind) = par.A1*log10(d_3d(ind)) + ...
                    par.B1 + ...
                    par.C1*log10(CenterFrequency) + ...
                    par.D1*log10(bp(ind).^2 + (hBS - hMS(ind)).^2);
                sf_sigma_tmp(ind) = par.sig1;
            end
            
            ind = ~ind;
            if sum(ind) > 0
                loss_1(ind) = par.A2*log10(d_3d(ind)) + ...
                    par.B2 + ...
                    par.C2*log10(CenterFrequency) + ...
                    par.D2*log10(bp(ind).^2 + (hBS - hMS(ind)).^2);
                sf_sigma_tmp(ind) = par.sig2;
            end
            
            loss_2 = par.A3*log10(d_3d) + ...
                par.B3 + ...
                par.C3*log10(CenterFrequency) + ...
                par.D3*(hMS - 1.5);
            
            use_loss_1 = bsxfun(@gt, loss_1, loss_2);
            loss(use_loss_1) = loss_1(use_loss_1);
            sf_sigma(use_loss_1) = sf_sigma_tmp(use_loss_1);
            
            use_loss_2 = ~use_loss_1;
            loss(use_loss_2) = loss_2(use_loss_2);
            sf_sigma(use_loss_2) = par.sig3;
            
        case 'urban_los'
            % From WINNER+ D5.3 pp. 74
            
            hBS = 1.5;
            hMS = 1.5;
            
            % Calculate the breakpoint
%             G = par.B1 + par.C1*log10(CenterFrequency) + par.D1*log10(hBS) + par.E1*log10(mean(hMS));
%             H = par.B2 + par.C2*log10(CenterFrequency) + par.D2*log10(hBS) + par.E2*log10(mean(hMS));
%             bp  = 10^( (H-G)/( par.A1-par.A2 ) );
            bp = 4*(hBS-1)*(hMS-1)*CenterFrequency*1e9/3e8;

            ind = d_2d<=bp;
            if sum(ind)>0
                loss(ind) = par.A1*log10(d_2d(ind)) + par.B1 + par.C1*log10(CenterFrequency)...
                    + par.D1*log10(hBS) + par.E1*log10(hMS(ind)) + par.F1*hMS(ind);
                sf_sigma(ind) = par.sig1;
            end
            
            ind = ~ind;
            if sum(ind)>0
                loss(ind) = par.A2*log10(d_2d(ind)) + par.B2 + par.C2*log10(CenterFrequency)...
                    + par.D2*log10(hBS) + par.E2*log10(hMS(ind)) + par.F2*hMS(ind);
                sf_sigma(ind) = par.sig2;
            end
            
        case 'freeway_los'
            % From WINNER+ D5.3 pp. 74
            
            hBS = 1.5;
            hMS = 1.5;
            
            % Calculate the breakpoint
%             G = par.B1 + par.C1*log10(CenterFrequency) + par.D1*log10(hBS) + par.E1*log10(mean(hMS));
%             H = par.B2 + par.C2*log10(CenterFrequency) + par.D2*log10(hBS) + par.E2*log10(mean(hMS));
%             bp  = 10^( (H-G)/( par.A1-par.A2 ) );
            bp = 4*(hBS-1)*(hMS-1)*CenterFrequency*1e9/3e8;

            ind = d_2d<=bp;
            if sum(ind)>0
                loss(ind) = par.A1*log10(d_2d(ind)) + par.B1 + par.C1*log10(CenterFrequency)...
                    + par.D1*log10(hBS) + par.E1*log10(hMS(ind)) + par.F1*hMS(ind);
                sf_sigma(ind) = par.sig1;
            end
            
            ind = ~ind;
            if sum(ind)>0
                loss(ind) = par.A2*log10(d_2d(ind)) + par.B2 + par.C2*log10(CenterFrequency)...
                    + par.D2*log10(hBS) + par.E2*log10(hMS(ind)) + par.F2*hMS(ind);
                sf_sigma(ind) = par.sig2;
            end
            
            
        case 'urban_nlos'
            % From WINNER+ D5.3 pp. 74
            
            hBS = 1.5;
            hMS = 1.5;
            
            % Calculate the breakpoint
            
%             G = par.B1 + par.C1*log10(CenterFrequency) + par.D1*log10(hBS) + par.E1*log10(mean(hMS));
%             H = par.B2 + par.C2*log10(CenterFrequency) + par.D2*log10(hBS) + par.E2*log10(mean(hMS));
%             bp  = 10^( (H-G)/( par.A1-par.A2 ) );
            bp = 4*(hBS-1)*(hMS-1)*CenterFrequency*1e9/3e8;
            
            % Calculate the PL_LOS(d_1)
            ind = d1_2d<=bp;
            if sum(ind)>0
                PL_LOS_d1(ind) = par.A1*log10(d1_2d(ind)) + par.B1 + par.C1*log10(CenterFrequency)...
                    + par.D1*log10(hBS) + par.E1*log10(hMS(ind)) + par.F1*hMS(ind);
            end
            
            ind = ~ind;
            if sum(ind)>0
                PL_LOS_d1(ind) = par.A2*log10(d1_2d(ind)) + par.B2 + par.C2*log10(CenterFrequency)...
                    + par.D2*log10(hBS) + par.E2*log10(hMS(ind)) + par.F2*hMS(ind);
            end
            
            % Calculate the PL_LOS(d_2)
            ind = d2_2d<=bp;
            if sum(ind)>0
                PL_LOS_d2(ind) = par.A1*log10(d2_2d(ind)) + par.B1 + par.C1*log10(CenterFrequency)...
                    + par.D1*log10(hBS) + par.E1*log10(hMS(ind)) + par.F1*hMS(ind);
            end
            
            ind = ~ind;
            if sum(ind)>0
                PL_LOS_d2(ind) = par.A2*log10(d2_2d(ind)) + par.B2 + par.C2*log10(CenterFrequency)...
                    + par.D2*log10(hBS) + par.E2*log10(hMS(ind)) + par.F2*hMS(ind);
            end
            
            % Calculate the PL(d1,d2)
            n_j = max((2.8-0.0024*d1_2d),1.84);
            PL_d1_d2 = par.A3*PL_LOS_d1 + par.B3 + par.C3*n_j + par.D3*n_j*log10(d2_2d) + par.E3*log10(CenterFrequency)+par.F3;
            
            % Calculate the PL(d2,d1)
            n_j = max((2.8-0.0024*d2_2d),1.84);
            PL_d2_d1 = par.A3*PL_LOS_d2 + par.B3 + par.C3*n_j + par.D3*n_j*log10(d1_2d) + par.E3*log10(CenterFrequency)+par.F3;
            
            % Calculate the PL
            loss = min(PL_d1_d2, PL_d2_d1);
            sf_sigma(:) = par.sig3;



            
        otherwise
            tmp = ones(size(sf_sigma));
            sf_sigma = h_parset.scenpar.SF_sigma + ...
                h_parset.scenpar.SF_delta * ...
                log10( h_parset.simpar.center_frequency / 1e9 );
            sf_sigma = tmp * sf_sigma;
    end
else
    tmp = ones(size(sf_sigma));
    sf_sigma = h_parset.scenpar.SF_sigma + ...
        h_parset.scenpar.SF_delta * ...
        log10( h_parset.simpar.center_frequency / 1e9 );
    sf_sigma = tmp * sf_sigma;
end

loss = reshape(loss,n_snapshots,n_rx,n_tx);

if use_track
    sf_sigma = mean(sf_sigma);
end

% The shadow fading might change with distance. Hence, if
% the value did change, we have to rescale the values from
% the map.
SF_sigma_scenpar = h_parset.scenpar.SF_sigma;
if SF_sigma_scenpar == 0
    scale_sf = 1;
else
    scale_sf = sf_sigma / SF_sigma_scenpar;
end

end
