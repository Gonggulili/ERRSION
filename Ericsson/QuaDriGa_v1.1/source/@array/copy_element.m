function copy_element( h_array, source, target )
%COPY_ELEMENT Creates a copy of an antenna element
%
%   COPY_ELEMENT(source,target) copies the element from source to target.
%   Source must be an index of the array object. The value must be scalar,
%   integer and greater than 0 and it can not exceed the array size.
%   Target can be any number > 0.
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input arguments
if ~( size(source,1) == 1 && isnumeric(source) && all(size(source) == [1 1]) ...
        &&  all( mod(source,1)==0 ) && min(source) > 0 && max(source)<=h_array.no_elements )
    error('??? "source" must be scalar, integer > 0 and can not exceed array size')
end

if ~( size(target,1) == 1 && isnumeric(target) ...
        &&  all( mod(target,1)==0 ) && min(target) > 0)
    error('??? "target" must be integer > 0')
end

if target > h_array.no_elements
    h_array.no_elements = max(target);
end

% Copy the data from source to target
target = setdiff(target,source);
for n=1:numel(target)
    h_array.element_position(:,target(n)) = h_array.element_position(:,source);
    h_array.Fa(:,:,target(n)) = h_array.Fa(:,:,source);
    h_array.Fb(:,:,target(n)) = h_array.Fb(:,:,source);
    
    if numel( h_array.Fc ) ~= 1
        h_array.Fc(:,:,target(n)) = h_array.Fc(:,:,source);
    end
end

end

