function han = visualize( obj, tx , rx, show_names , create_new_figure )
%VISUALIZE Plot the layout
%
%   VISUALIZE creates a plot of the layout.
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input arguments
if exist('tx','var') && ~isempty( tx )
    if isempty(tx)
        tx = 1:obj.no_tx;
    elseif ~( size(tx,1) == 1 && isnumeric(tx) ...
            &&  all( mod(tx,1)==0 ) && min(tx) > 0 && max(tx)<=obj.no_tx )
        error('??? "tx" must be integer > 0 and can not exceed array size')
    end
else
    tx = 1:obj.no_tx;
end

if exist('rx','var') && ~isempty( rx )
    if isempty(rx)
        rx = 1:obj.no_rx;
    elseif ~( size(rx,1) == 1 && isnumeric(rx) ...
            &&  all( mod(rx,1)==0 ) && min(rx) > 0 && max(rx)<=obj.no_rx )
        error('??? "rx" must be integer > 0 and can not exceed array size')
    end
else
    rx = 1:obj.no_rx;
end

if exist('show_names','var') && ~isempty( show_names )
    if ~( all(size(show_names) == [1,1]) && isnumeric(show_names) )
        error('??? "show_names" must be scalar and numeric.')
    end
else
    show_names = 1;
end

if exist('create_new_figure','var') && ~isempty( create_new_figure )
    if ~( all(size(create_new_figure) == [1,1]) && isnumeric(create_new_figure) )
        error('??? "create_new_figure" must be scalar and numeric.')
    end
else
    create_new_figure = true;
end

% Get the ground directions
for n=1:obj.no_rx
    if isempty( obj.track(n).ground_direction )
        obj.track(n).compute_directions;
    end
end

if create_new_figure
    han = figure('Position',[ 100 , 100 , 1000 , 700]);
end

plot3( obj.tx_position(1,tx(1)),obj.tx_position(2,tx(1)),obj.tx_position(3,tx(1)),...
    '+r','Linewidth',3,'Markersize',15 );
hold on

pos_tx_ant = zeros(3,obj.tx_array(tx(1)).no_elements);
for n=1:3
    pos_tx_ant(n,:) = obj.tx_position(n,tx(1)) + obj.tx_array(tx(1)).element_position(n,:);
end
plot3( pos_tx_ant(1,:),pos_tx_ant(2,:),pos_tx_ant(3,:),...
    '^r','Linewidth',1,'Markersize',10 );

seg_pos = obj.track(rx(1)).initial_position;
plot3( seg_pos(1,:),seg_pos(2,:),seg_pos(3,:),...
    'ob','Linewidth',2,'Markersize',10 );

pos_rx_ant = obj.rx_array(rx(1)).element_position;
deg = obj.track(rx(1)).ground_direction(1);

Rz = makehgtform('zrotate', deg);
pos_rx_ant = Rz(1:3, 1:3) * pos_rx_ant;
for n=1:3
    pos_rx_ant(n,:) = obj.track(rx(1)).initial_position(n,1) + pos_rx_ant(n,:);
end

plot3( pos_rx_ant(1,:),pos_rx_ant(2,:),pos_rx_ant(3,:),...
    'vb','Linewidth',1,'Markersize',10 );

pos = obj.track(rx(1)).positions;
for n=1:3
    pos(n,:) = pos(n,:) + obj.track(rx(1)).initial_position(n);
end
plot3( pos(1,:),pos(2,:),pos(3,:) );


for m = 1:numel(tx)
    plot3( [obj.tx_position(1,tx(m)),obj.tx_position(1,tx(m))] , ...
        [obj.tx_position(2,tx(m)),obj.tx_position(2,tx(m))],...
        [0,obj.tx_position(3,tx(m))] ,'--r' );
    
    if m>1
        plot3( obj.tx_position(1,tx(m)),obj.tx_position(2,tx(m)),obj.tx_position(3,tx(m)),...
            '+r','Linewidth',3,'Markersize',15 );
        
        pos_tx_ant = zeros(3,obj.tx_array(tx(m)).no_elements);
        for n=1:3
            pos_tx_ant(n,:) = obj.tx_position(n,tx(m)) + obj.tx_array(tx(m)).element_position(n,:);
        end
        plot3( pos_tx_ant(1,:),pos_tx_ant(2,:),pos_tx_ant(3,:),...
            '^r','Linewidth',1,'Markersize',10 );
    end
    
    if show_names
        text( obj.tx_position(1,tx(m)), obj.tx_position(2,tx(m)), obj.tx_position(3,tx(m)),...
            obj.tx_name{tx(m)}) ;
    end
end


for m = 1:numel(rx)
    si = obj.track(rx(m)).segment_index;
    pos = obj.track(rx(m)).positions;
    for n=1:3
        pos(n,:) = pos(n,:) + obj.track(rx(m)).initial_position(n);
    end
    sii = floor( si + 0.3*([si(2:end),obj.track(rx(m)).no_snapshots ]-si ) );
    
    if m>1
        pos_rx_ant = obj.rx_array(rx(m)).element_position;
        deg = obj.track(rx(m)).ground_direction(1);
        
        Rz = makehgtform('zrotate', deg);
        pos_rx_ant = Rz(1:3, 1:3) * pos_rx_ant;
        for n=1:3
            pos_rx_ant(n,:) = obj.track(rx(m)).initial_position(n,1) + pos_rx_ant(n,:);
        end
        plot3( pos(1,:),pos(2,:),pos(3,:) );
        plot3( pos_rx_ant(1,:),pos_rx_ant(2,:),pos_rx_ant(3,:),...
            'vb','Linewidth',1,'Markersize',10 );
    end
    plot3( pos(1,si),pos(2,si),pos(3,si) ,'ob','Linewidth',2,'Markersize',10 );
    
    if show_names == 1
        for o = 1:numel(sii)
            tmp = obj.track(rx(m)).scenario(:,o);
            tmp = regexprep(tmp,'_','\\_');
            text(pos(1,sii(o)),pos(2,sii(o)),pos(3,sii(o)),tmp) ;
        end
    elseif show_names == 2
        text(obj.track(rx(m)).initial_position(1),...
            obj.track(rx(m)).initial_position(2),...
            obj.track(rx(m)).initial_position(3),...
            obj.track(rx(m)).name) ;
    end
end

hold off

grid on
box on
view(0, 90);

xlabel('x-coord in [m]');
ylabel('y-coord in [m]');
zlabel('z-coord in [m]');

legend('Tx-Position','Tx-Antenna','Rx-Position','Rx-Antenna','Rx-Track',...
    'Location','NorthEastOutside')

a = axis;
A = a(2)-a(1);
B = a(4)-a(3);
if A<B
    a(1) = a(1) - 0.5*(B-A);
    a(2) = a(2) + 0.5*(B-A);
else
    a(3) = a(3) - 0.5*(A-B);
    a(4) = a(4) + 0.5*(A-B);
end

a(1) = min( a(1) , -1 );
a(2) = max( a(2) , 1 );
a(3) = min( a(3) , -1 );
a(4) = max( a(4) , 1 );

a = a*1.1;
axis(a);


end



