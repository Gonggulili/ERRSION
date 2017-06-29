function [ len, dist ] = get_length( obj )
%GET_LENGTH Calculates the length of the track in [m]
%
%   len = GET_LENGTH calculate the length of a track in m.
%
%   Output variables:
%   	"len":
%           Length of a track in [m]
%
%   	"dist":
%           Distance of each position (snapshot) from the start of the
%           track in [m] 
%
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

no_trk = numel( obj );

len = zeros( 1,no_trk );
dist = cell( 1,no_trk );

for m = 1:no_trk
    p = obj(m).positions;
    for n=1:3
        p(n,:) = p(n,:) - p(n,1);
    end
    
    if isempty( obj(m).Plength ) || nargout == 2
        
        dist{m} = zeros(1,obj(m).no_snapshots);
        for n=2:obj(m).no_snapshots
            dist{m}(n) = sqrt(sum(( p(:,n) - p(:,n-1) ).^2));
        end
        dist{m} = cumsum(dist{m});
        len(m) = dist{m}(end);
        obj(m).Plength = len(m);
        
    else
        len(m) = obj(m).Plength;
    end
end

if no_trk == 1
    dist = dist{1};
else
   len = reshape( len , size(obj) );
   dist = reshape( dist , size(obj) );
end

end

