function uncompress( h_array )
%UNCOMPRESS Uncompresses an array
%
%   See: array.compress
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if h_array.iscompressed
    Fa = h_array.Fa;
    Fb = h_array.Fb;
    
    h_array.Pind = [];
    
    h_array.PFa = Fa;
    h_array.PFb = Fb;
    
    h_array.precision = 'double';
end

end
