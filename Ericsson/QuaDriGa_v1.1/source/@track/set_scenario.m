function set_scenario( h_track, scenario, probability,...
    seg_length_min, seg_length_mu, seg_length_std )
%SET_SCENARIO Assigns random scenarios and creates segments
%
% SET_SCENARIO( scenario, probability, seg_length_min, seg_length_mu, seg_length_std  ) 
%   This function can be used to create segments along the trajectory and
%   assign scenarios to the segments. 
%
%   Input variables:
%   	"scenario":
%           a cell array of scenario-names. Each scenario (synonym for
%           propagation environment) is described by a string (e.g.
%           "MIMOSA_16-25_LOS" or "WINNER_SMa_C1_NLOS"). A list of
%           supported scenarios can be obtained by calling
%           "parameter_set.supported_scenarios". The scenario parameters
%           are stored in the configuration folder "config" in the QuaDRiGa
%           main folder. The filenames (e.g. "MIMOSA_16-25_LOS.conf") also
%           serves as scenario name.       
%
%   	"probability":
%           the probability for which the scenario occurs. This parameter
%           must be a vector of the same length as there are scenarios.
%           Probabilities must be specified in between 0 and 1. The sum of
%           the probabilities must be 1. By default (or when "probability"
%           is set to "[]"), each scenario is equally likely.     
%
%   	"seg_length_min":
%           the minimal segment length in [m]. The default is 10 m.
%
%   	"seg_length_mu":
%           the median segment length in [m]. The default is 30 m. 
%
%   	"seg_length_std":
%           the standard deviation of the street length in [m]. The default
%           is 12 m.
%
%   If there are less than 3 input arguments (i.e. only "scenario" and/or
%   "probability" is given), then no segments will be created. To create
%   segments with the default settings call "SET_SCENARIO(scenario,[],[])".
%
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if ~exist( 'scenario' , 'var' )
    error('"scenario" is not given.')
end
if ischar(scenario)
    scenario_list_old = scenario;
    scenario = cell(1,1);
    scenario{1} = scenario_list_old;
end


% Parse input arguments
supported_scenarios = parameter_set.supported_scenarios;

% Add file names to list of supported scenarios
for i_scen = 1:numel( scenario )
    if numel(scenario{i_scen}) > 5 && ...
            ~isempty( regexp( scenario{i_scen}(end-4:end) , '.conf', 'once' ) ) && ...
            regexp( scenario{i_scen}(end-4:end) , '.conf' ) == 1
        file = dir(scenario{i_scen});
        if ~isempty(file)
            supported_scenarios{end+1} = scenario{i_scen};
        end
    end
end

str = ' not found. Supported scenarios are: ';
no = numel(supported_scenarios);
for n = 1:no
    str = [str,supported_scenarios{n}];
    if n<no
        str = [str,', '];
    end
end
no_scenario_list = numel( scenario );
for n = 1 : no_scenario_list
    if ~any( strcmpi(scenario{n},supported_scenarios) )
        error(['??? Scenario type "',scenario{n},'"',str]);
    end
end


% Parse "probability"
if ~exist( 'probability' , 'var' ) || isempty(probability)
    probability = ones(no_scenario_list,1)./no_scenario_list;
elseif numel(probability) ~= numel(scenario)
    error('??? number of elements in "probability" must match the number of elements in "scenario".');
elseif ~( isnumeric(probability) && isreal(probability) &&...
        all( probability >= 0) && all( probability <= 1) )
    error('??? "probability"has wrong format.');
end
probability = probability(:) ./ sum( probability(:) );
cum_probability = cumsum( probability );

create_segments = false;
if nargin > 3
    create_segments = true;
    
    % Parse "seg_length_min"
    if ~exist( 'seg_length_min' , 'var' ) || isempty( seg_length_min )
        seg_length_min = 10;
    elseif ~(isnumeric( seg_length_min ) && seg_length_min>=0 &&...
            isreal( seg_length_min ) && all(size(seg_length_min) == [1 1]))
        error('??? "seg_length_min"  has wrong format');
    end
    
    % Parse "seg_length_mu"
    if ~exist( 'seg_length_mu' , 'var' ) || isempty( seg_length_mu )
        seg_length_mu = 30;
    elseif ~(isnumeric( seg_length_mu ) && seg_length_mu>=0 &&...
            isreal( seg_length_mu ) && all(size(seg_length_mu) == [1 1]))
        error('??? "seg_length_mu"  has wrong format');
    end
    
    % Parse "seg_length_std"
    if ~exist( 'seg_length_std' , 'var' ) || isempty( seg_length_std )
        seg_length_std = 12;
    elseif ~(isnumeric( seg_length_std ) && seg_length_std>=0 &&...
            isreal( seg_length_std ) && all(size(seg_length_std) == [1 1]))
        error('??? "seg_length_mu"  has wrong format');
    end
end


for i_track = 1:numel( h_track )
    
    % There are several scenarios, one for each tx in the layout
    % Determine the number of Tx:
    n_tx = size( h_track(i_track).scenario ,1);
    
    if create_segments
                
        [ trk_length , dist ] = h_track(i_track).get_length;
        no_snapshots = h_track(i_track).no_snapshots;
        
        pos = 0;
        no_segments = 1;
        while pos < trk_length
            
            % Get the length of the current segment
            seg_length = randn*seg_length_std + seg_length_mu;
            while seg_length < seg_length_min
                seg_length = randn*seg_length_std + seg_length_mu;
            end
            ind = find( pos + seg_length < dist,1 );
            pos = dist( ind );
            
            if ind < no_snapshots     % Exception for the last segment
                no_segments = h_track(i_track).no_segments + 1;
                h_track(i_track).no_segments = no_segments;
                h_track(i_track).segment_index( no_segments ) = ind;
                
                if ~isempty( scenario )
                    for i_tx = 1:n_tx
                        i_scen = find( rand < cum_probability,1 );
                        h_track(i_track).scenario{ i_tx,no_segments-1 } = scenario{i_scen};
                    end
                end
            end
        end
        % Set last scenario
        if ~isempty( scenario )
            for i_tx = 1:n_tx
                i_scen = find( rand < cum_probability,1 );
                h_track(i_track).scenario{ i_tx,no_segments } = scenario{i_scen};
            end
        end
        
    else
        % If segments are already given, assign random scenarios
        if ~isempty( scenario )
            for i_segment = 1:h_track(i_track).no_segments
                for i_tx = 1:n_tx
                    i_scen = find( rand < cum_probability,1 );
                    h_track(i_track).scenario{ i_tx,i_segment } = scenario{i_scen};
                end
            end
        end
    end
end

end
