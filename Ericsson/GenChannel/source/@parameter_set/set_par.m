function set_par( h_parset, name, value )
%SET_PAR Sets the parameters of all objects in 'parameter_set' arrays.
%
% SET_PAR( name,value ) sets all values of the parameter specified by
% 'name'  of the 'parameter-set'-array to the given value.
% Example: set_par( 'ds' , 1e-9 ) sets all ds-values to one ns.
%
% Input:
%   "name":
%       The field name that should be altered. Supported are:
%       'scenpar', 'plpar',
%       'ds', 'kf', 'sf', 'asD', 'asA', 'esD', 'esA',
%       'LSP_xcorr_matrix'
%       'samples_per_meter', and 'map_extension'.
%
%   "value":
%       The value that should be assigned.
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


supported_types = { 'scenpar','plpar','ds','kf','sf','asD','asA','esD','esA',...
    'samples_per_meter','map_extension','LSP_xcorr_matrix'};

if ~( ischar(name) && any( strcmpi(name,supported_types)) )
    str = 'Name not found; supported types are: ';
    no = numel(supported_types);
    for n = 1:no
        str = [str,supported_types{n}];
        if n<no
            str = [str,', '];
        end
    end
    error(str);
end


switch name
    case {'scenpar','plpar'}
        
        for i_parset = 1 : numel(h_parset)
            if i_parset == 1
                h_parset(i_parset).(name) = value;
            else
                h_parset(i_parset).(['P',name]) = value;
            end
        end
        
    case 'LSP_xcorr_matrix'
        
        if ~all( size(value) == [7,7] )
            error('QuaDRiGa:Parameter_set:wrongInputValue',...
                'LSP_xcorr_matrix must be of size [7x7]')
        end
        
        if ~all( diag(value) == ones(7,1) )
            error('QuaDRiGa:Parameter_set:wrongInputValue',...
                'LSP_xcorr_matrix must have all ones on main diagonal')
        end
        
        if any(any( value>1 | value<-1 ))
            error('QuaDRiGa:Parameter_set:wrongInputValue',...
                'Entries in LSP_xcorr_matrix be in between -1 and 1')
        end
        
        for i_parset = 1 : numel(h_parset)
            sp = h_parset(i_parset).scenpar;
            
            sp.ds_kf = value( 1,2 ) ;
            sp.ds_sf = value( 1,3 ) ;
            sp.asD_ds = value( 1,4 ) ;
            sp.asA_ds = value( 1,5 ) ;
            sp.esD_ds = value( 1,6 ) ;
            sp.esA_ds = value( 1,7 ) ;
            sp.sf_kf = value( 2,3 ) ;
            sp.asD_kf = value( 2,4 ) ;
            sp.asA_kf = value( 2,5 ) ;
            sp.esD_kf = value( 2,6 ) ;
            sp.esA_kf = value( 2,7 ) ;
            sp.asD_sf = value( 3,4 ) ;
            sp.asA_sf = value( 3,5 ) ;
            sp.esD_sf = value( 3,6 ) ;
            sp.esA_sf = value( 3,7 ) ;
            sp.asD_asA = value( 4,5 ) ;
            sp.esD_asD = value( 4,6 ) ;
            sp.esA_asD = value( 4,7 ) ;
            sp.esD_asA = value( 5,6 ) ;
            sp.esA_asA = value( 5,7 ) ;
            sp.esD_esA = value( 6,7 ) ;
            
            h_parset(i_parset).scenpar = sp;
        end
        
    otherwise
        
        if ~( all(size(value) == [1,1]) && isnumeric(value) && isreal(value) )
            error('??? "value" must be scalar and numeric')
        end
        
        for n = 1:numel(h_parset)
            eval( [ 'h_parset(n).',name,' = ones( 1,numel(h_parset(n).',name,') )*value;' ] );
        end
        
end

end

