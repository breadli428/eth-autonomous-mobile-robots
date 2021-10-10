function [nodes, parents, solution] = rrt(map, goalIdx, parameters, startIdx)
% Rapidly exploring random trees for finding a path between the goal and start points in map.
% Shorthands
L_ = parameters.max_branch_length; % default 5
P_ = parameters.goal_fixation_probability; % default 0.05
D_ = parameters.debug_plot; % default false
R_ = 0; % goal radius

% the tree is a list of nodes, and corresponding ids of each node's parent
% as the tree grows, nodes should have dimension   (n_nodes, 2)
%                    parents should have dimension (n_nodes, 1)
nodes = startIdx;
parents = 1; % the first nodes' parent is itself
paths = {[]}; % paths between nodes

% turn the map image (y,x) into an (x,y) occupancy grid
isfree = map';
map_size = size(isfree);

% Initialize plotting
if D_
    hold off
    imshow(isfree')
    set(gca,'YDir','normal')
    hold on
    scatter(startIdx(1), startIdx(2), 'MarkerFaceColor',[0 .7 .7])
    scatter(goalIdx(1), goalIdx(2), 'MarkerFaceColor',[.7 0 0])
end

while 1
    % Randomly sample a point 'qrand' in (x,y) in [(1, 1), (map_size_x, map_size_y)]
    % With probability P_, select the goal coordinates as qrand.
    %STARTRM
    if (rand(1) < P_)
    %ENDRM
    %STARTUNCOMMENT
%     if (TODO)
    %ENDUNCOMMENT
        qrand = goalIdx;
    else
        %STARTRM
        % Generates two integers in [(1, 1), (map_size_x, map_size_y)]
        qrand = round(rand(1,2) .* (map_size-1))+1;
        %ENDRM
        %STARTUNCOMMENT
%         qrand = TODO
        %ENDUNCOMMENT
    end
    % Find closest node 'qnear' to sampled point 'qrand'
    %STARTRM
    [~, qnear_id] = min(vecnorm(qrand - nodes, 2, 2));
    %ENDRM
    %STARTUNCOMMENT
%     qnear_id = TODO % the index in nodes of the nearest node.
    %ENDUNCOMMENT
    qnear = nodes(qnear_id,:);
    if all(qnear == qrand)
        continue
    end
    % Grow branch according to local rules 
    points = path_(qnear, qrand, isfree, L_);
    if size(points,1) == 1
        continue
    end
    qnew = points(end,:);
    % Add the new node to the tree
    nodes = [nodes; qnew];
    parents = [parents; qnear_id];
    paths = [paths; points];
    % Debug plot
    if D_
        plot(points(:,1)', points(:,2)')
        scatter(qnew(1), qnew(2))
        drawnow
    end
    % Success condition
    if norm(qnew - goalIdx) <= R_
        break
    end
end

% Extract the solution from the tree
shortest_path_coarse = nodes(end,:);
shortest_path_fine = flipud(paths{end});
parent = parents(end);
while 1
    shortest_path_coarse = [shortest_path_coarse; nodes(parent,:)];
    path_to_previous_node = flipud(paths{parent});
    path_to_previous_node = path_to_previous_node(2:end,:);
    shortest_path_fine = [shortest_path_fine; path_to_previous_node];
    parent = parents(parent);
    if parent == 1
        shortest_path_coarse = [shortest_path_coarse; nodes(1,:)];
        path_to_previous_node = flipud(paths{parent});
        path_to_previous_node = path_to_previous_node(2:end,:);
        shortest_path_fine = [shortest_path_fine; path_to_previous_node];
        break
    end
end
shortest_path_coarse = flipud(shortest_path_coarse);
shortest_path_fine = flipud(shortest_path_fine);

% Plot solution
if D_
    for i = 1:size(shortest_path_coarse,1)
        a = shortest_path_coarse(i,:);
        scatter(a(1), a(2), 'MarkerEdgeColor',[.7 0 0])
    end
    plot(shortest_path_fine(:,1)', shortest_path_fine(:,2)', 'r')
end

% Return solution
solution = shortest_path_fine;
end

function [points] = path_(start, goal, isfree, max_length)
    % Finds the shortest path between two pixels, where only
    % up/down/left/right movement is allowed.
    % The result is an array of points (n_points, n_dimensions)
    % containing the starting point and all subsequent points in the path.
    %STARTRM
    eye_ = eye(2); % utility array
    
    candidate = start;
    points = candidate;
    % Grow path between start and goal. As long as a step towards the goal
    % is possible, take that step. The step taking us closer to the goal is
    % preferred.
    for i = 1:max_length
        d = goal - candidate;
        if all(d == 0)
            return
        end
        [~, strong_dim] = max(abs(d));
        weak_dim = 2 + 1 - strong_dim;
        strong_step = eye_(strong_dim,:) .*  sign(d);
        weak_step = eye_(weak_dim,:) .* sign(d);
        
        step = candidate + strong_step;
        if isfree(step(1), step(2))
            candidate = step;
            points = [points; candidate];
            continue
        end
        
        if norm(weak_step) == 0
            break
        end
        
        step = candidate + weak_step;
        if isfree(step(1), step(2))
            candidate = step;
            points = [points; candidate];
            continue
        end
        break
    end
    %ENDRM
    %STARTUNCOMMENT
%     points = [start; TODO]
    %ENDUNCOMMENT
end
