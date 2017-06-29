function correct_overlap( h_track , overlap )
%CORRECT_OVERLAP Correct position of the segment start to account for overlap
%
%   After the channel coefficients are calculated, adjacent segments can be
%   merged into one time-continuous output. The merger assumes that the
%   merging interval happens at the end of one segment, before a new
%   segments starts. In a reality,  however, the scenario change happens in
%   the middle of the overlapping part (and not at the end of it).      
%
%   This function corrects the position of the segment start to account for
%   that. The parameter "overlap" indicates, how long the overlapping part
%   of the segment is. It can have values in between 0 and 1.    
% 
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input arguments
check = true;
if nargin < 2
    overlap = 0.5;
    check = false;
end

if check
    if ~( isnumeric(overlap) && all(size(overlap) == [1 1]) && isreal(overlap) ...
            && overlap<=1 && overlap>=0 )
        error('??? Overlap must be scalar, and in between 0 and 1')
    end
end

if numel(h_track) > 1
    % Do for each element in array
    for n=1:numel(h_track)
        h_track(n).correct_overlap( overlap );
    end
else
    % Only do if there are more than on segment
    if h_track.no_segments > 1
        % Get a list of the segment indices
        seg_ind_old = [ h_track.segment_index , h_track.no_snapshots ];
        seg_ind_new = zeros( size(  h_track.segment_index ));
        seg_ind_new(1) = 1;

        for n = 1:h_track.no_segments-1
            seg_length = seg_ind_old(n+1) - seg_ind_new(n);
            overlapping = seg_length * overlap;
            new_seg_length = round( seg_length + 0.588235*overlapping );
            
            if seg_ind_new(n) + new_seg_length < seg_ind_old(n+2)
                seg_ind_new(n+1) = seg_ind_new(n) + new_seg_length;
            else
                seg_ind_new(n+1) = seg_ind_old(n+1);
            end
        end
        
        % Assign new segment index
        h_track.segment_index = seg_ind_new;
                
%         overlap_point = zeros( size(  h_track.segment_index ));
%         overlap_point(1) = 1;
%         for n = 1:h_track.no_segments-1
%             seg_length = seg_ind_new(n+1) - seg_ind_new(n);
%             overlapping = seg_length * overlap;
%             overlap_point(n+1) = round( seg_ind_new(n) + seg_length - 0.5*overlapping );
%         end
                
    end
end

end