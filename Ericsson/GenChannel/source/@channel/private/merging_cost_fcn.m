function [ mse , weight ] = merging_cost_fcn( oA , oB , p , d , ds_target , ramp  )
%MERGING_COST_FCN Calculates costs for the channel merger
%
%   When merging the channel coefficients of adjacent segments in a track,
%   we do not want the delay spread to change. The distribution of the DS
%   is given in the parameter tables (prameter_set.scenpar) as log-normal
%   distributed. However, during merging, new taps ramp up and old taps
%   ramp down. This will of course have an effect on the delay spread. Each
%   subsegment of an overlapping path will thus have a different delay
%   spread depending on which taps are ramped up and down. "oA" and "oB"
%   contain a permutation of taps. The first one is always the LOS. The
%   later ones are for NLOS. Ramping is done in the order given in "oA" and
%   "oB". I.e. at first, tab number two [oA(2)] ramps down and at the same
%   time [oB(2)] ramps up. Then [oA(3)] ramps down and at the same time
%   [oB(3)] ramps up. All taps not having a partner will be ramped up/down
%   without a counterpart.             
%
%   This function calculates the delay spread for each subsegment of the
%   merging interval and returns MSE compared to a linear ramp between the
%   DS of the first segment and the second.  
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


% Determine the lengths of the vectors
L1 = numel(oA);
L2 = numel(oB);

if L1 >= L2
    no_subseg = L2-1;
else
    no_subseg = L1-1;
end

% Here, we calculate the weight matrix for each supsegment.
% This will later be used to estimate the DS during merging.

weight = zeros(no_subseg , L1+L2);
tmp = eye(no_subseg) * 0.5;
one = ones(no_subseg);
ind = 2:no_subseg+1;
weight( : , oA( ind )    ) = triu( one ) - tmp;
weight( : , oB( ind )+L1 ) = tril( one ) - tmp;

if L1>L2
    weight( : , oA( no_subseg+2:L1 ) ) = 1-ramp(:,ones(1,L1-no_subseg-1));
elseif L2>L1
    weight( : , oB( no_subseg+2:L2 )+L1 ) = ramp(:,ones(1,L2-no_subseg-1));
end

weight( :,1 ) = 1-ramp;
weight( :,L1+1 ) = ramp;

% We calculate the DS for each subsegment
p = weight .* p( ones(1,no_subseg) , : );
ds = sqrt( p*d.^2 - (p*d).^2 );

% The MSE for the DS
mse = sqrt( sum( ( ds_target - ds ).^2 ) );
