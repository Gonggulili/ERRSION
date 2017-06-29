function subtracks = get_subtrack( obj , i_segment )
%GET_SUBTRACK Splits the track in subtracks for each segment (state)
%
%   subtracks = GET_SUBTRACK( i_segment ) returns the subtracks for the given
%   segment indices. When no input argument is provided, all subtracks are
%   returned.
% 
%   After defining segments along the track, one needs the subtrack
%   that corresponds only to one segment to perform the channel calculation.
%   This new track can consist of two segments. The first segment contains
%   the positions from the previous segment, the second from the current.
%   This is needed to generate overlapping channel segments for the merging
%   process.
%
%   Input and output variables:
%   	"i_segment":
%           A list of indices indicating which subtracks should be
%           returned. By default, all subtracks are returned.
%
%   	"subtracks":
%           A vector of track objects corresponding to the number of
%           segments. 
%
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if nargin == 1
   i_segment = 1:obj.no_segments;
   check = false;
else
    check = true;
end

% Parse the input variables
if check
    if ~( any(size(i_segment) == 1) && isnumeric(i_segment) ...
            && isreal(i_segment) && all(mod(i_segment,1)==0) && all(i_segment > 0) )
        error('??? "i_segment" must be integer and > 0')
    elseif max(i_segment) > obj.no_segments
        error('??? "i_segment" exceeds number of entries in object')
    end
end

% Make a copy of the parameters
par = obj.par;
names = {'ds','kf','pg','asD','asA','esD','esA','xpr'};

for seg=1:numel(i_segment)
    segment = i_segment(seg);
    
    % Select the part of the track that corresponds to the given segment
    if obj.no_segments == 1
        % The current track has only one segment
        if obj.closed == 1
            ind = 1:obj.no_snapshots-1;
        else
            ind = 1:obj.no_snapshots;
        end
        seg_ind = 1;
        scen_ind = 1;
        
    elseif segment == 1
        % The current track has more than one segment, and the returned
        % segment is the first one.
        if obj.closed == 1
            ind = [ obj.segment_index( obj.no_segments ) : obj.no_snapshots-1 ,...
                1:obj.segment_index( segment+1 )-1 ];
            seg_ind = [ 1 , obj.no_snapshots-obj.segment_index( obj.no_segments )+1 ];
            scen_ind = [ obj.no_segments , 1 ];
        else
            ind = 1 : obj.segment_index( segment+1 )-1;
            seg_ind = 1;
            scen_ind = 1;
        end
        
    elseif segment == obj.no_segments
        % The current track has more than one segment, and the returned
        % segment is the last one.
        if obj.closed == 1
            ind = obj.segment_index( segment-1 ) : obj.no_snapshots-1;
        else
            ind = obj.segment_index( segment-1 ) : obj.no_snapshots;
        end
        seg_ind = [ 1 , obj.segment_index( segment ) - obj.segment_index( segment-1 ) + 1 ];
        scen_ind = [ segment-1 , segment ];
        
    else
        % The current track has more than one segment, and the returned
        % segment neither the first, nor the last one.
        ind = obj.segment_index( segment-1 ) : obj.segment_index( segment+1 )-1;
        seg_ind = [ 1 , obj.segment_index( segment ) - obj.segment_index( segment-1 ) + 1 ];
        scen_ind = [ segment-1 , segment ];
        
    end
    
    % Create new track with the corresponding data
    tr = track;
    tr.name                     = [obj.name,'_',num2str(segment)];
    tr.positions                = obj.positions( :,ind );
    tr.segment_index            = seg_ind;
    if tr.no_segments == 2
        sp = tr.positions( :, seg_ind(2) );
        tr.scenario = obj.scenario(:,scen_ind([2,2]));
    else
        sp = tr.positions( :, 1 );
        tr.scenario = obj.scenario(:,scen_ind);
    end
    
    if ~isempty( obj.ground_direction )
        tr.ground_direction         = obj.ground_direction(ind);
        tr.height_direction         = obj.height_direction(ind);
    end
    
    if ~isempty(par)
        % Init struct
        out_par = struct('ds',[],'kf',[],'pg',[],'asD',[],'asA',[],'esD',[],'esA',[],'xpr',[]);
        
        % Copy the data
        for n = 1:8
            val = par.(names{n});
            if ~isempty( val )
                if n==2 || n==3
                    out_par.( names{n} ) = val(:,ind);
                else
                    out_par.( names{n} ) = val(:,scen_ind);
                end
            end
        end
        
        % Save to subtrack
        tr.par_nocheck = out_par;
    end
    
    % Set the initial position
    tr.initial_position = obj.initial_position + sp;
    for n=1:3
        tr.positions(n,:) = tr.positions(n,:) - sp(n);
    end
    
    % Append to subtracks-list
    subtracks(seg) = tr;
end
end

