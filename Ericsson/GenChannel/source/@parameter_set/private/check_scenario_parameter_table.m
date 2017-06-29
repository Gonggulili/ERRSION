function [out, ok] = check_scenario_parameter_table(value, check)
%CHECK_SCENARIO_PARAMETER_TABLE Checks the scenario table
%
%   CHECK_SCENARIO_PARAMETER_TABLE is a private function of the parameter_set
%   class. It can not be called from outside the class members and is thus not
%   directly accessible to the user. However, when the scenpar structure of a
%   parameter_set object is edited, check_scenario_parameter_table is immediately
%   called to check the validity of the new parameters. Changing the Winner table
%   (scenpar) in a loop can thus be very slow. This functions also calculates the
%   cross-correlation matrix (LSP_xcorr_matrix) and checks if it is positive
%   definite (LSP_matrix_isOK). If this matrix is not positive definite, then
%   the generation of correlated LSPs will fail (update_parameters will return an
%   error).
%
% QuaDRiGa Copyright (C) 2011-2016 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if nargin == 1
    check = true;
end

if check
    
    t_int = {'NumClusters'};
    
    t_positive = {'r_DS','DS_lambda','AS_D_lambda','AS_A_lambda',...
        'ES_D_lambda','ES_A_lambda','SF_lambda','KF_lambda',};
   
    t_positive_or_zero = {'PerClusterAS_D','PerClusterAS_A','PerClusterES_D',...
        'PerClusterES_A','LNS_ksi','LOS_scatter_radius'};
    
    t_real = {'xpr_mu','xpr_sigma','xpr_gamma','xpr_delta',...
        'DS_mu','DS_sigma','DS_gamma','DS_delta',...
        'AS_D_mu','AS_D_sigma','AS_D_gamma','AS_D_delta',...
        'AS_A_mu','AS_A_sigma','AS_A_gamma','AS_A_delta',...
        'ES_D_mu','ES_D_mu_A','ES_D_mu_min','ES_D_sigma','ES_D_gamma','ES_D_delta',...
        'ES_A_mu','ES_A_sigma','ES_A_gamma','ES_A_delta',...
        'SF_sigma','SF_delta',...
        'KF_mu','KF_sigma','KF_gamma','KF_delta',...
        'NumClusters_gamma'};
    
    t_abs1 = { 'asD_ds','asA_ds','asA_sf','asD_sf','ds_sf','asD_asA','asD_kf',...
        'asA_kf','ds_kf','sf_kf','esD_ds', 'esA_ds','esA_sf','esD_sf','esD_esA',...
        'esD_asD','esD_asA','esA_asD','esA_asA','esD_kf','esA_kf'};
    
    names = fieldnames(value);
    if numel(names) ~= 69
        error('??? Wrong number of fields in "scenpar".');
    end
    for n = 1:numel(names)
        v = value.(names{n});
        
        if ismember( names(n) , t_int ) &&...
                ~( all(size(v) == [1 1]) && isnumeric(v) && isreal(v) && mod(v,1)==0 && v > 0 )
            error(['??? "',names{n},'" must be integer, scalar and > 0'])
        end
        
        if ismember( names(n) , t_positive ) &&...
                ~( all(size(v) == [1 1]) && isnumeric(v) && isreal(v) && v > 0 )
            error(['??? "',names{n},'" must be real, scalar and > 0'])
        end
        
        if ismember( names(n) , t_positive_or_zero ) &&...
                ~( all(size(v) == [1 1]) && isnumeric(v) && isreal(v) && v >= 0 )
            error(['??? "',names{n},'" must be real, scalar and >= 0'])
        end
        
        if ismember( names(n) , t_real ) &&...
                ~( all(size(v) == [1 1]) && isnumeric(v) && isreal(v) )
            error(['??? "',names{n},'" must be real and scalar'])
        end
        
        if ismember( names(n) , t_abs1 ) &&...
                ~( all(size(v) == [1 1]) && isnumeric(v) && isreal(v) && v<=1 && v>=-1 )
            error(['??? "',names{n},'" must be real, scalar and have values from -1 to 1'])
        end
        
    end
    
end

% In order to create correlated LSPs, we need to create a correlation
% matrix. THis is done here, when the input of the parameter table is
% OK. The output is saved as a property.

a = value.ds_kf;          % delay spread vs k-factor
b = value.ds_sf;          % delay spread vs shadowing std
c = value.asD_ds;         % departure AS vs delay spread
d = value.asA_ds;         % arrival AS vs delay spread
e = value.esD_ds;         % departure ES vs delay spread
f = value.esA_ds;         % arrival ES vs delay spread
g = value.sf_kf;          % shadowing std vs k-factor
h = value.asD_kf;         % departure AS vs k-factor
k = value.asA_kf;         % arrival AS vs k-factor
l = value.esD_kf;         % departure DS vs k-factor
m = value.esA_kf;         % arrival DS vs k-factor
n = value.asD_sf;         % departure AS vs shadowing std
o = value.asA_sf;         % arrival AS vs shadowing std
p = value.esD_sf;         % departure ES vs shadowing std
q = value.esA_sf;         % arrival ES vs shadowing std
r = value.asD_asA;        % departure AS vs arrival AS
s = value.esD_asD;        % departure ES vs departure AS
t = value.esA_asD;        % arrival ES vs departure AS
u = value.esD_asA;        % departure ES vs arrival AS
v = value.esA_asA;        % arrival ES vs arrival AS
w = value.esD_esA;        % departure ES vs arrival ES


% Cross-correlation matrix
% Order of rows and columns is ds, kf, sf, asD, asA, esD, esA
out = [ 1  a  b  c  d  e  f ;...
    a  1  g  h  k  l  m ;...
    b  g  1  n  o  p  q ;...
    c  h  n  1  r  s  t ;...
    d  k  o  r  1  u  v ;...
    e  l  p  s  u  1  w ;...
    f  m  q  t  v  w  1 ];

ok = true;

if check
    % check for positive definiteness
    [~, p] = chol(out, 'lower');
    if p > 0
        warning('MATLAB:parameter_set',...
            ['LSP_xcorr_matrix is not positive-definite. \n',...
            'Generation of correlated LSPs will fail. \n',...
             'Eigenvalues:      ',sprintf('%1.2f  ',eig(out))]);
        ok = false;
    end
end
