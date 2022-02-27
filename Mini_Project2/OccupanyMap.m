clear all; close all;

%room configuration
roomspec;

%generate occupancy map
roomMap = roomOccupancyGrid(10,10, colobj);

%show occupancy map
figure(1)
show(roomMap);
view(-90,90);
axis([-1 11 -1 11]);
grid;



roomshow(colobj, 2);