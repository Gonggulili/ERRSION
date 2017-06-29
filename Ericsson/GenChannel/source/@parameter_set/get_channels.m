function [ h_channel, h_cb  ] = get_channels( h_parset )
%GET_CHANNELS Generates channel coefficients
%
%   c = GET_CHANNELS generates the channel coefficients. This is
%   the main function of the channel builder. 
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


[N,O] = size(h_parset);

m = 1;
for n=1:N
    for o=1:O
        if h_parset(n,o).no_positions > 0
            h_cb(m) = channel_builder( h_parset(n,o) );
            m = m+1;
        end
    end
end

h_channel = h_cb.get_channels;

end

