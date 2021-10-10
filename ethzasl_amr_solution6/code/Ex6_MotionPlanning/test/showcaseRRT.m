close all
fileDir = fileparts(mfilename('fullpath'));

map = loadMapFromImage( [fileDir, '/../maps/simple_100x100.png'] );
map = flipud(1 + map);

start = [65, 82];
goal = [80, 15];

parameters.max_branch_length = 5;
parameters.debug_plot = true;
parameters.goal_fixation_probability = 0.05;

[~,~, solution] = rrt(map, goal, parameters, start);
%STARTRM

disp("Press any key for next example");
pause;

close all
map = loadMapFromImage( [fileDir, '/../maps/simple_100x100.png'] );
map = flipud(1 + map);

start = [65, 82];
goal = [80, 15];

parameters.max_branch_length = 12;
parameters.debug_plot = true;
parameters.goal_fixation_probability = 0.05;

[~,~, solution] = dynamic_rrt(map, goal, parameters, start);

disp('Press any key for next example');
pause;

close all
map = loadMapFromImage( [fileDir, '/../maps/random_map.png'] );
map = flipud(1 + map(:,:,1));

start = [30, 41];
goal = [274, 264];

parameters.max_branch_length = 15;

[~,~, solution] = rrt(map, goal, parameters, start);
%ENDRM
