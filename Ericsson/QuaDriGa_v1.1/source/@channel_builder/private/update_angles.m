function [ ang_upd, spread_upd ] = update_angles( spread, ang, pow, fixed, range )
% RESCALE_ANGLES Iterativley updates the angles to match a given angular spread
%
% Input variables:
%   spread      Vector of angular spreads [ N x 1 ]
%   ang         Matrix of angles in rad [ N x L ]
%   pow         Matrix of linear power values [ N x L ]
%   fixed       True / False matrix indicating which angles are fixed [ N x L ]
%   range       Allowed rangle for the angles (min, max) [ 2 x 1 ]
%
% Output variables:
%   ang_upd     Matrix of updated angles in rad [ N x L ]
%   spread_upd  Vector of updated angular spreads [ N x 1 ]
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

N    = numel(spread);
L    = size(ang,2);
oL   = ones(1,L);

% Calculate the average angle for the fixed angles
avg_ang_fixed = zeros( N,1 );
for n = 1:N
    ii = fixed(n,:);
    avg_ang_fixed(n) = angle( sum( pow(n,ii).*exp( 1j*ang(n,ii) ) , 2 ) ); % [rad]
end

% Index list of the angles that need to be chaecked  / updated
upd = true(N,1);
spread_upd = NaN( N,1 );

cnt_outer_loop = 0;
while any( upd ) && cnt_outer_loop < 20
    cnt_outer_loop = cnt_outer_loop + 1;
    
    % Check if there are angles outside specified range
    replace_ang = ang < range(1) | ang > range(2);
    replace_ang( fixed ) = false;
    
    cnt_inner_loop = 0;
    while any( replace_ang( : ) ) && cnt_inner_loop < 20
        cnt_inner_loop = cnt_inner_loop + 1;
        
        % Generate new angles spread around the fixed angle range
        nnn = find( any( replace_ang,2 ) );
        for nn = 1:numel(nnn)
            n = nnn(nn);
            no_new_ang = sum( replace_ang(n,:) );
            new_ang = spread(n) * randn( 1,no_new_ang ) + avg_ang_fixed(n);
            ang( n,replace_ang(n,:) ) = new_ang;
        end
        
        replace_ang = ang < range(1) | ang > range(2);
        replace_ang( fixed ) = false;
    end
    
    % Calculate the angular spreads
    spread_upd( upd ) = qf.calc_angular_spreads( ang(upd,:), pow(upd,:) );
    
    % Target accuracy is 0.001 rad
    upd = abs( spread-spread_upd ) > 1e-3;
    
    scale_new   = spread ./ spread_upd;
    scale_new( ~upd ) = 1;
    scale       = ones( N,1 );
    step_size   = scale_new - scale;
    
    chn = upd;
    
    cnt_inner_loop = 1;
    while any( chn ) && cnt_inner_loop < 20
        cnt_inner_loop = cnt_inner_loop + 1;
        
        scale_new( chn ) = scale( chn ) + step_size( chn );
        
        % Caculate the scaled angles
        ang_new = ang( chn,: ) .* scale_new( chn,oL );
        ang_new( fixed( chn,: ) ) = ang( fixed & chn( :,oL ) );
        
        % Update the angular spread
        spread_new = qf.calc_angular_spreads( ang_new, pow(chn,:) );
        
        % Target accuracy is 0.001 rad
        ok = abs( spread( chn ) - spread_new ) < 1e-3 | abs(step_size( chn ) ) < 1e-4;
        
        % Values with increases match
        bt = abs( spread( chn ) - spread_new ) < abs( spread( chn ) - spread_upd( chn ) );
        
        chn_ind = find( chn );
        scale( chn_ind( bt|ok ) )   = scale_new( chn_ind( bt|ok ) );
        spread_upd( chn_ind( bt|ok ) ) = spread_new( bt|ok );
        
        % Reduce step size and change direction
        step_size( chn_ind(~bt) ) = -0.37 * step_size( chn_ind(~bt) );
        step_size( chn_ind(ok) )  = 0;
        
        ang(  chn_ind(ok),: ) = ang_new( ok,: );
        
        chn( chn_ind(ok) ) = false;
    end
    
    if any( chn )
        ang( chn,: ) = ang_new( ~ok,: );
        spread_upd( chn ) = spread_new( ~ok );
    end
    
end

ang_upd = ang;

end
