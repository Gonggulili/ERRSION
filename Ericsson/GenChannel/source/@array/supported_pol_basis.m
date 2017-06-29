function out = supported_pol_basis( value )
% SUPPORTED_POL_BASIS Provides a list of the currently supported polarization bases
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

supported_types = {'cartesian','az-el','el-az','polar-spheric'};
if nargin == 1
    if ischar(value) && any( strcmpi(value,supported_types))
        out = true;
    else
        out = false;
    end
else
    out = supported_types;
end

end
