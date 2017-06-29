function [ scenarios, file_names, file_dir ] = supported_scenarios( parse_shortnames )
%SUPPORTED_SCENARIOS Returns a list of supported scenarios
%
%   The optional parameter "parse_shortnames" can disable the shortname
%   parsing. This is significantly faster.
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input variables.
if nargin > 0
    parse_shortnames = logical(parse_shortnames(1));
else
    parse_shortnames = true;
end

% Parse config files in current folder
files_tmp = dir(['*.conf']);

files = [];
m = 1;
for n = 1:numel(files_tmp)
    files(m).name = files_tmp(n).name;
    files(m).dir  = '';
    m = m + 1;
end

% Parse files in config-folder
information = what('@parameter_set');
config_folder = regexprep('/home/ubuntu/gongyi/Ericsson/GenChannel/source/@parameter_set','@parameter_set',['config',filesep]);
files_tmp = dir([config_folder,'*.conf']);

for n = 1:numel(files_tmp)
    files(m).name = files_tmp(n).name;
    files(m).dir  = config_folder;
    m = m + 1;
end

nfiles = numel(files);

scenarios = cell(1,1000);
if parse_shortnames
    file_names = cell(1,1000);
    file_dir   = cell(1,1000);
end

m = 0;
for n = 1 : nfiles
    m = m + 1;
    [~, scenarios{m}] = fileparts( files(n).name );
    file_names{m} = files(n).name;
    file_dir{m}   = files(n).dir;
    
    if parse_shortnames
        % Read the file content to find short names for the scenario
        file = fopen( [ files(n).dir , files(n).name ] ,'r');
        lin = fgets(file);
        while ischar(lin)
            p1 = regexp(lin,'[ \t]*ShortName[ \t]*=[ \t]*[A-Za-z0-9_]+');
            if p1 == 1
                short_name = ...
                    regexp( lin,...
                    '(?<=^[ \t]*ShortName[ \t]*=[ \t]*)[A-Za-z0-9_]+',...
                    'match');
                
                if any( strcmp( short_name{1} , scenarios(1:m) ) )
                    error( ['ShortName "',short_name{1},'" in "',files(n).name,'" is ambiguous.'] )
                end
                
                m = m+1;
                scenarios{m}  = short_name{1};
                file_names{m} = files(n).name;
                file_dir{m}   = files(n).dir;
            end
            lin = fgets(file);
        end
        fclose(file);
    end
end

scenarios = scenarios(1:m);
if parse_shortnames
    file_names = file_names(1:m);
    file_dir = file_dir(1:m);
end

end

