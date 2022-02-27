%choose your algorithms here

algor = 4;
% 1 -> RRT
% 2 -> PRM
% 3 -> Potential Field
% 4 -> bug algorithm


if algor == 1
    RRT;
elseif algor == 2
    PRM;
elseif algor == 3
    Potential_field;
elseif algor == 4
    bug_algorithm;
end