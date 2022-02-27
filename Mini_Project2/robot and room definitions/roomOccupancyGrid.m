%OccupancyGrid
function map=roomOccupancyGrid(width, height, colobj)
    res = 1/10;

    roomMap = binaryOccupancyMap(width, height, 10);

    setOccupancy(roomMap, [0 0], ones(1,width*10));
    setOccupancy(roomMap, [width-res 0], ones(height*10,1));
    setOccupancy(roomMap, [0 0], ones(height*10,1));
    setOccupancy(roomMap, [0 height-res], ones(1,width*10));
    
    for i=5:length(colobj.obj)
        x = colobj.pos{i}(1);
        y = colobj.pos{i}(2);

        if(colobj.type(i) == 1)%if this is a box type object
            bx = x - (colobj.size{i}(1)/2);
            by = y - (colobj.size{i}(2)/2);
            w = colobj.size{i}(2)*10;
            h = colobj.size{i}(1)*10;
            obj = ones(w, h);
            setOccupancy(roomMap, [bx by], obj);
        %if this is a cylinder or sphere type object
        elseif(colobj.type(i) == 2 || colobj.type(i) == 3)
            r = colobj.size{i}(1);
            theta = linspace(0,2*pi);
            xc = x+r*cos(theta);
            yc = y+r*sin(theta);
            setOccupancy(roomMap, [xc' yc'], 1);
                
        end
    
    end
    map = roomMap;
end
