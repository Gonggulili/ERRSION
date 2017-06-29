function out = copy_objects(obj)
%COPY_OBJECTS is a modified version of the standard physical copy function
%
%   While the standard copy command creates new physical objects for each
%   element of obj (in case obj is an array of object handles), copy_objects checks
%   whether there are object handles pointing to the same object and keeps
%   this information.
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

for n = 1:numel(obj)
    ind = obj(n).eq(obj);
    ind(n) = false;
    ind = find(ind,1);
    
    if ind < n
        out(n) = out(ind);
    else
        out(n) = obj(n).copy;
    end
end

end