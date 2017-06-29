function [ el, az ] = maxi( P, elevation_grid, azimuth_grid, res )

% Interpolation targets
res = 1/ceil(1/res);
u   = 0:res:1;
nu  = numel(u);

nEl = numel(elevation_grid);
nAz = numel(azimuth_grid);

if abs(azimuth_grid(1)+pi) < 1e-6 && abs(azimuth_grid(end)-pi) < 1e-6
    wrap_around = true;
else
    wrap_around = false;
end

% Get initial Values
[~,ind] = max( P(:) );
[iEl,iAz] = ind2sub(size(P),ind);

% Get the neighboring elevation positions
pEl1 = int16( [ iEl-2 ; iEl-1 ; iEl ; iEl+1 ] );
pEl2 = int16( [ iEl-1 ; iEl ; iEl+1 ; iEl+2 ] );
if pEl1(1) < 1
    pEl2 = pEl2-pEl1(1)+1;
    pEl1 = pEl1-pEl1(1)+1;
elseif pEl2(end) > nEl
    pEl1 = pEl1 - (pEl2(end)-nEl);
    pEl2 = pEl2 - (pEl2(end)-nEl);
end
pEl1( pEl1 > nEl ) = nEl;
pEl2( pEl2 > nEl ) = nEl;
pEl1( pEl1 < 1 )   = 1;
pEl2( pEl2 < 1 )   = 1;
vEl = elevation_grid(pEl1(2:4));

% Get the neighboring azimuth positions
pAz1 = int16( [ iAz-2 ; iAz-1 ; iAz ; iAz+1 ] );
pAz2 = int16( [ iAz-1 ; iAz ; iAz+1 ; iAz+2 ] );
if pAz1(1) < 1
    if wrap_around
        pAz1(pAz1==0)  = nAz-1;
        pAz1(pAz1==-1) = nAz-2;
        pAz2(pAz2==0)  = nAz-1;
    else
        pAz2 = pAz2-pAz1(1)+1;
        pAz1 = pAz1-pAz1(1)+1;
    end
elseif pAz2(end) > nAz
    if wrap_around
        pAz1(pAz1==nAz+1) = 2;
        pAz2(pAz2==nAz+1) = 2;
        pAz2(pAz2==nAz+2) = 3;
    else
        pAz1 = pAz1 - (pAz2(end)-nAz);
        pAz2 = pAz2 - (pAz2(end)-nAz);
    end
end
pAz1( pAz1 > nAz ) = nAz;
pAz2( pAz2 > nAz ) = nAz;
pAz1( pAz1 < 1 )   = 1;
pAz2( pAz2 < 1 )   = 1;
vAz = unwrap(azimuth_grid( pAz1(2:4) ));

% Kernel for the spline interpolation function
K = 0.5 .* ...
    [0,  2,  0,  0; ...
    -1,  0,  1,  0; ...
     2, -5,  4, -1; ...
    -1,  3, -3,  1].';

% Get the offset between the actual user position and the next point
% on the map.
uK = K * [ ones(1,numel(u)) ; u ; u.^2 ; u.^3];

% Interpolate
Pi = zeros(2*nu-1);
Pi( 1:nu , 1:nu  ) = uK.' * P(pEl1,pAz1) * uK;
Pi( 1:nu ,nu:end ) = uK.' * P(pEl1,pAz2) * uK;
Pi(nu:end, 1:nu  ) = uK.' * P(pEl2,pAz1) * uK;
Pi(nu:end,nu:end ) = uK.' * P(pEl2,pAz2) * uK;

% Find the maximum
[~,ind] = max( Pi(:) );
[iEli,iAzi] = ind2sub(size(Pi),ind);

% Determine final angle
if vEl(1) ~= vEl(3)
    v  = vEl(1) : mean(diff(vEl))/(nu-1) : vEl(3);
    el = v(iEli);
else
    el = vEl(2);
end

if vAz(1) ~= vAz(3)
    v = vAz(1) : mean(diff(vAz))/(nu-1) : vAz(3);
    az = angle(exp(1j*v(iAzi)));
else
    az = vAz(2);
end

end