function [ aod,eod,aoa,eoa,delay ] = get_subpath_angles( h_cb,i_mobile )
%GET_SUBPATH_ANGLES Generate subpaths and perform random coupling (private)
%
%   GET_SUBPATH_ANGLES generates the 20 subpaths around the each path and
%   randomly couples the subpaths on the Tx- and Rx side. 
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

use_ground_reflection = h_cb.par.simpar.use_ground_reflection;

n_clusters = h_cb.NumClusters;                     % no. clusters
cpl = h_cb.subpath_coupling(:,:,:,i_mobile);       % Subpath couping

% The 20 offset angles for the sub-paths
offset = [0.0447 0.1413 0.2492 0.3715 0.5129 0.6797 0.8844 1.1481 1.5195 2.1551];
offset = [ offset , -offset ]*0.017453292519943;  % in rad

% The per cluster angular spread scaling coefficients
c_aod = h_cb.par.scenpar.PerClusterAS_D;    % in deg
c_aoa = h_cb.par.scenpar.PerClusterAS_A;    % in deg
c_eod = h_cb.par.scenpar.PerClusterES_D;    % in deg
c_eoa = h_cb.par.scenpar.PerClusterES_A;    % in deg

% Reserve some memory for the output
aod = zeros( 1,n_clusters,20 );    
aoa = aod;
eod = aod;
eoa = aod;

% We get the subpath angles at the Rx position
p = zeros(20,4);
for i_cluster=1:n_clusters
    p(:,1) = c_aod * offset( cpl(:,1,i_cluster) );
    p(:,2) = c_aoa * offset( cpl(:,2,i_cluster) );
    p(:,3) = c_eod * offset( cpl(:,3,i_cluster) );
    p(:,4) = c_eoa * offset( cpl(:,4,i_cluster) );
    
    if i_cluster == 1 || ( i_cluster == 2 && use_ground_reflection )
        p = zeros(20,4);
    end

    aod(1,i_cluster,:) = h_cb.AoD( i_mobile,i_cluster ) + p(:,1);
    aoa(1,i_cluster,:) = h_cb.AoA( i_mobile,i_cluster ) + p(:,2);
    eod(1,i_cluster,:) = h_cb.EoD( i_mobile,i_cluster ) + p(:,3);
    eoa(1,i_cluster,:) = h_cb.EoA( i_mobile,i_cluster ) + p(:,4);
end

% Calculate delays
if nargout == 5
    n_snapshots = h_cb.par.rx_track(i_mobile).no_snapshots;
    if h_cb.par.simpar.use_absolute_delays
        r_0 = h_cb.par.rx_track(i_mobile).initial_position - h_cb.par.tx_position;
        D = sqrt(sum(r_0.^2)) / h_cb.par.simpar.speed_of_light;
        delay = repmat( h_cb.taus(i_mobile,:)+D , n_snapshots , 1  );
    else
        delay = repmat( h_cb.taus(i_mobile,:) , n_snapshots , 1  );
    end
end

end
