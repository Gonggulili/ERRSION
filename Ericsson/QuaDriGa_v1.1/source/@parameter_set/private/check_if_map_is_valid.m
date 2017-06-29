function valid = check_if_map_is_valid( scenpar_old , scenpar )
%CHECK_IF_MAP_IS_VALID Checks it the map is still valid after parameter change
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


names = fieldnames( scenpar_old );

% If any of the following parameters are changes, no new maps need to be
% generated. 
uncritical = {'NumClusters','r_DS','PerClusterAS_D','PerClusterAS_A','PerClusterES_D',...
    'PerClusterES_A','LOS_scatter_radius',...
    'LNS_ksi','xpr_mu','xpr_sigma'};

critical = ~ismember(names,uncritical);

O = cell2mat( struct2cell(scenpar_old) );
N = cell2mat( struct2cell(scenpar) );

if any( O-N ~= 0 & critical )
    valid = false;
else
    valid = true;
end

