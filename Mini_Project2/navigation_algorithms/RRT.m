clear all; close all;

roomspec;

map = roomOccupancyGrid(10,10, colobj);

ss = stateSpaceSE2;
sv = validatorOccupancyMap(ss);

sv.Map = map;

sv.ValidationDistance = 0.5;
ss.StateBounds = [map.XWorldLimits; map.YLocalLimits; [-pi pi]];
resolution = 10;

planner = plannerRRT(ss, sv);
planner.MaxConnectionDistance = 0.5;

start = [1 1 0];
goal = [8 8 0];

rng(100,"twister");

[pthObj,solnInfo] = plan(planner,start,goal);

figure(1)
show(map);
hold on
plot(solnInfo.TreeData(:,1),solnInfo.TreeData(:,2),'.-'); % tree expansion
plot(pthObj.States(:,1),pthObj.States(:,2),'r-','LineWidth',2) % draw path
view([-90 90]);