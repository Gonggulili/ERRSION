function generate_initial_xpr( h_cb )
%GENERATE_INITIAL_XPR Generates the initial XPR values

L = h_cb.NumClusters;                      % no. taps
N = h_cb.par.no_positions;                 % no. positions
oL = ones(1,L-1);
oP = ones(1,20);

log10_f_GHz = log10( h_cb.par.simpar.center_frequency / 1e9 );

mu     = h_cb.par.scenpar.xpr_mu;
gamma  = h_cb.par.scenpar.xpr_gamma;
sigma  = h_cb.par.scenpar.xpr_sigma;
delta  = h_cb.par.scenpar.xpr_delta;

mu = mu + gamma * log10_f_GHz;
sigma = sigma + delta * log10_f_GHz;
sigma( sigma<0 ) = 0;

if h_cb.par.simpar.use_polarization_rotation == 0 ||...
        h_cb.par.simpar.use_polarization_rotation == 3
    
    % Use polarization model from WINNER
    xpr_subpath = randn(N,L-1,20)*sigma + mu;
    h_cb.xpr    = cat( 2, ones(N,1,20) * Inf , xpr_subpath );
    
else
    %% Linear Polarization
    
    % We get the mean value from the parameter set and add an additional
    % spread for the NLOS clusters. This spread is equal to the original
    % XPR-sigma. 
    xpr_mu      = 10*log10( h_cb.par.xpr.' );
    xpr_subpath = randn(N,L-1,20)*sigma + xpr_mu(:,oL,oP);
    h_cb.xpr    = cat( 2, ones(N,1,20) * Inf , xpr_subpath );
    
end

%% Circular Polarization

% Set the XPR for the circular polarization
xpr_mu      = 10*log10( h_cb.par.xpr.' );
xpr_cluster = randn(N,L-1)*sigma + xpr_mu * oL ;
xpr_cluster = 10.^( 0.1*xpr_cluster );

rand_sign   = 2*( randi(2,size(xpr_cluster))-1.5 );

kappa_cluster = [zeros( N,1 ) , rand_sign.*acot( sqrt(xpr_cluster) ) ];
h_cb.kappa  = kappa_cluster(:,:,oP);


%% Random Polarization

h_cb.random_pol = (rand(4,L*20,N)*2*pi - pi);

end
