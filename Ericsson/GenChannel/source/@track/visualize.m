function visualize( obj )
%VISUALIZE Plots the track
%
%   VISUALIZE creates a plot of the track.
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

p3 = [];
for m = 1:numel(obj)    
    p3 = [p3 ,  unique(obj(m).positions(3,:))+obj(m).initial_position(3)];
end
if numel(unique(p3))==1
    space = false;
else
    space = true;
end

figure('Position',[ 100 , 100 , 1000 , 700]);

for m = 1:numel(obj)
    si = obj(m).segment_index;
    pos = obj(m).positions;
    for n=1:3
        pos(n,:) = pos(n,:) + obj(m).initial_position(n);
    end
    sii = floor( si + 0.3*([si(2:end),obj(m).no_snapshots ]-si ) );
    
    if space
        plot3( pos(1,:),pos(2,:),pos(3,:) );
        if m==1
            hold on
        end
        plot3( pos(1,si),pos(2,si),pos(3,si) ,'or' );
        
        tmp = obj(m).scenario;
        tmp = regexprep(tmp,'_','\\_');
        text(pos(1,sii),pos(2,sii),pos(3,sii),tmp) ;
    else
        plot( pos(1,:),pos(2,:) );
        if m==1
            hold on
        end
        plot( pos(1,si),pos(2,si),'or' );
        for n = 1:obj(m).no_segments
            tmp = obj(m).scenario(:,n);
            tmp = regexprep(tmp,'_','\\_');
            text(pos(1,sii(n)),pos(2,sii(n)),tmp) ;
        end
    end
    
end
hold off

grid on
title('Track layout')
xlabel( 'x-coord in [m]' )
ylabel( 'y-coord in [m]' )
zlabel( 'z-coord in [m]' )
legend('Track','Segment start','Location','NorthEastOutside')

a = axis;
A = a(2)-a(1);
B = a(4)-a(3);
if A<B
    a(1) = a(1) - 0.5*(B-1);
    a(2) = a(2) + 0.5*(B-A);
else
    a(3) = a(3) - 0.5*(A-B);
    a(4) = a(4) + 0.5*(A-B);
end

a(1) = min( a(1) , -1 );
a(2) = max( a(2) , 1 );
a(3) = min( a(3) , -1 );
a(4) = max( a(4) , 1 );

if space
    a(5) = min( a(5) , -1 );
    a(6) = max( a(6) , 1 );
end
a = a*1.1;
axis(a);



end



