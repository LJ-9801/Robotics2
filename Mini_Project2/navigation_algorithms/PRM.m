clear all; close all;

roomspec;

map = roomOccupancyGrid(10,10, colobj);



prm = mobileRobotPRM(map, 5);

startloc = [1 1];
endloc = [9 5];
path = findpath(prm, startloc, endloc);
show(prm);
view([-90 90]);