function h_parset = generate_map_winner( h_parset, delta, sigma, vb_dots )
% GENERATE_MAP Generates the correlation map for the LSPs
%
% This function generates the maps for the correlates large scale
% parameters by digital filtering a random normal distributed process.
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if ~exist('vb_dots','var')
    vb_dots = 0;
    verbose = false;
elseif vb_dots == 0
    verbose = false;
else
    verbose = true;
end

% Get map dimensions
no_map_y = h_parset.map_size(2);
no_map_x = h_parset.map_size(1);

% Initialize output map
no_delta = numel(delta);

% Correct the input
delta = 0.5*delta;

% Repeat for each delta-value
m0=0;
for n = 1 : no_delta
    if verbose; m1=ceil(n/no_delta*vb_dots); if m1>m0; for m2=1:m1-m0; fprintf('o'); end; m0=m1; end; end;
    
    if sigma(n) == 0
        % If there is no variance in the parameter, we don't need to create
        % any autocorrelation map. This saves time and memory.
        mapi = 0;
        
    else
        
        % The offset
        do = ceil( delta(n)*4 );
        
        % Limit the maximum filter length to 1000 to stop MATLAB from freezing
        do = min( [do,1000] );
        
        % The filter function for the elements in x/y-direction
        d = 0:do;
        h = exp(-1 * d./delta( ones(1,numel(d))*n  ));
        h = h ./ sqrt(delta(n));
        
        % The filter-function for the diagonal direction
        d2 = 0:sqrt(2):do;
        h2 = exp(-1 * d2./delta( ones(1,numel(d2))*n  ));
        h2 = h2 ./ sqrt(delta(n));
        
        % Extend the map
        nrow = no_map_y+2*do;
        ncol = no_map_x+do;
        
        rind = do+1 : do+no_map_y;
        cind = do+1 : do+no_map_x;
        
        % Initialize the map with random coefficients
        mapi = randn( nrow,ncol,'single' );
        
        % Calculate the x/y correlation
        mapi = filter( h,1,mapi,[],1 );
        mapi = filter( h,1,mapi,[],2 );
        
        % For debugging only
        %             imagesc(mapi)
        %             hold on
        %             plot(cind,ones(size(cind))*rind(end),'-k')
        %             plot(cind,ones(size(cind))*rind(1),'-k')
        %             plot([ones(size(rind))]*rind(1),rind,'-k')
        %             hold off
        %             pause
        
        % Extract the center part of the map and adjust the mu and sigma of the
        % normal distribution
        mapi = mapi( rind , cind );
        mu  = mean( reshape( mapi , 1 ,[] )  );
        sig = std( reshape( mapi , 1 ,[] )  );
        mapi = (mapi-mu)./sig;
        
    end
    
    switch n
        case 1
            h_parset.ds_map  = mapi;
        case 2
            h_parset.kf_map  = mapi;
        case 3
            h_parset.sf_map  = mapi;
        case 4
            h_parset.asD_map = mapi;
        case 5
            h_parset.asA_map = mapi;
        case 6
            h_parset.esD_map = mapi;
        case 7
            h_parset.esA_map = mapi;
    end
    
end

