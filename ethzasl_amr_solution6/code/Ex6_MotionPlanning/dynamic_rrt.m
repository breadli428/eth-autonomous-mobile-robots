function [nodes, parents, solution] = dynamic_rrt(map, goalIdx, parameters, startIdx)
% Function to perform RRT taking into account robot dynamics.
% This does not need to be filled in and will be shown in the solutions as
% a rough example for one of many possible extensions for the algorithm.
%STARTRM
L_ = parameters.max_branch_length; % default 5
A_ = 1; % length of extra collision checking (look ahead)
P_ = parameters.goal_fixation_probability; % default 0.05
D_ = parameters.debug_plot; % default false
N_ = 30; % number of trajectory samples
R_ = 3; % goal radius
C_ = 24; % angle quantization constant

% the tree is a list of nodes, and corresponding parents
nodes = [startIdx 0];
paths = {[]};
parents = 1;

% for some reasons images are indexed (y, x)
isfree = map';
map_size = size(isfree);

samples_nrot = createTrajectorySamples(N_, L_ * (1+A_), C_); 

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
    % Goal-biased random sampling
    if (rand(1) < P_)
        qrand = goalIdx(1:2);
    else
        qrand = round(rand(1,2) .* (map_size-1))+1;
    end
    d = nodes;
    d(:,1:2) = qrand - d(:,1:2); % x y differences
    d(:,3) = atan2(-d(:,1), d(:,2))- d(:,3); % angle difference
    d = mc_(d);
    % Find closest node to sampled point
    [~, qnear_id] = min(vecnorm(d, 2, 2));
    qnear = nodes(qnear_id,:);
    if all(qnear(1:2) == qrand)
        continue
    end
    qrand(3) = atan2(-d(qnear_id,1), d(qnear_id,2));
    % Extract samples with rotation corresponding to qnear
    angle_id = min(round(mod(qnear(3),(2*pi)) * C_ / (2*pi)) + 1, C_);
    samples = reshape(samples_nrot(angle_id,:,:,:),N_,[],3);
    % Grow branch according to local rules 
    points = path_(qnear, qrand, isfree, samples);
    if size(points,1) < 2
        continue
    end
    qnew = points(end,:);
    % Add new node
    if any(all(qnew == nodes,2))
        continue
    end
    nodes = [nodes; qnew];
    paths = [paths; points];
    parents = [parents; qnear_id];
    % Debug plot
    if D_
        plot(points(:,1)', points(:,2)')
        scatter(qnew(1), qnew(2))
        drawnow
    end
    % Success condition
    if norm(qnew(1:2) - goalIdx) <= R_
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

function [points] = path_(start, goal, isfree, samples)
    % sort trajectories in terms of middle of traj closest to goal
    D_ = size(samples,3); % dimensions
    
    samples_aligned = round(samples + reshape([start(1:2), 0],1,1,[]));
    [costs, closest] = min(vecnorm(mc_(samples_aligned(:,1:ceil(end/2),:)-reshape(goal,1,1,3)),2,3),[],2);
    [~,indices] = sort(costs,1);
    best_trajectory = [];
    for s = indices'
        optimal_length = closest(s) * 2;
        trajectory = reshape(samples_aligned(s,:,:),[],D_);
        if optimal_length < size(trajectory,1)
            trajectory = trajectory(1:optimal_length,:);
        end
        % check how far this trajectory can be taken
        prev_point = start;
        free_trajectory = [];
        for t = 1:size(trajectory,1)
            point = trajectory(t,:);
            if ~pixel_path_(prev_point(1:2), point(1:2), isfree) || t == size(trajectory,1) || all(point(1:2) == goal(1:2))
                if t > 1
                  free_trajectory = trajectory(1:ceil((t-1)/2),:);
                end
                break
            end
        end
        % compare it to the best estimate
        if isempty(best_trajectory)
            best_trajectory = free_trajectory;
            continue
        end
        if isempty(free_trajectory)
           break
        end
        best_cost = norm(mc_(best_trajectory(end,:)-goal));
        cost = norm(mc_(free_trajectory(end,:)-goal));
        % if no improvement give up
        if cost < best_cost
            best_trajectory = free_trajectory;
        else
            break
        end
    end
    
    points = [start; best_trajectory];
end

function [d_new] = mc_(d_)
   % Metric correction
   % d_ is a difference vector in x y theta space
   % or an array of vectors (last dimension is x y theta)
   W_ = [1 1 20]; % weight of each dimension in calculating node cost
   d_new = d_;
   if length(size(d_)) == 3
     d_new(:,:,3) = min(abs(d_new(:,:,3)), abs((2*pi) - abs(d_new(:,:,3))));
     d_new = d_new .* reshape(W_,1,1,3);
   else
     d_new(:,3) = min(abs(d_new(:,3)), abs((2*pi) - abs(d_new(:,3))));
     d_new = d_new .* W_;
   end
end

function [path_is_clear] = pixel_path_(start, goal, isfree)
    eye_ = eye(2); % utility array
    
    candidate = start;
    path_is_clear = 0;
    
    % Check that goal is in isfree
    if any(goal > size(isfree) | goal < 1)
        return
    end
    % Grow pixel-wise path between start and goal
    while 1
        d = goal - candidate;
        if all(d == [0 0])
            path_is_clear = 1;
            return
        end
        [~, strong_dim] = max(abs(d));
        weak_dim = 2 + 1 - strong_dim;
        strong_step = eye_(strong_dim,:) .*  sign(d);
        weak_step = eye_(weak_dim,:) .* sign(d);
        
        step = candidate + strong_step;
        if isfree(step(1), step(2))
            candidate = step;
            continue
        end
        
        if norm(weak_step) == 0
            break
        end
        
        step = candidate + weak_step;
        if isfree(step(1), step(2))
            candidate = step;
            continue
        end
        break
    end
    return
end

function [samples] = createTrajectorySamples(n_samples, trajectory_length, n_rotations)
    L_ = trajectory_length; % longest trajectory
    R_ = L_/pi; % tightest turn radius (should be >L_/pi to avoid >180deg turns)
    n_R_ = n_samples - mod(n_samples,2); % samples in radius dimension
    C_ = n_rotations;
    
    % Parametrize curve with constant velocity of 1 along perimeter
    dpdt = 1;
    p = dpdt:dpdt:L_;
    
    % Generate curves of constant r
    r_inv = linspace(-1/R_,1/R_,n_R_)';
    theta = p .* r_inv;
    xx = (cos(theta) - 1) ./ r_inv;
    yy = sin(theta) ./ r_inv;
    samples = cat(3,xx,yy,theta); % dimensions (i_traj, j_point, k_xytheta)
    %plot(samples(:,:,1)', samples(:,:,2)')
    
    % generate an array of rotated samples for each angle quantization.
    rot_samples = zeros(C_, n_R_, L_, 3);
    angles = linspace(0,(2*pi),C_+1);
    angles = angles(1:end-1); % ignore 2pi
    rot_samples(1,:,:,:) = reshape(samples,1,n_R_,L_,[]);
    for i = 2:length(angles)
        a = angles(i);
        R = [cos(a), -sin(a); sin(a), cos(a)];
        for j = 1:n_R_
            % rotate a single trajectory and add it to the collection.
            trajectory = reshape(samples(j,:,:),L_,[]);
            trajectory(:,1:2) = (R * trajectory(:,1:2)')';
            trajectory(:,3) = round(C_ / (2*pi) * (trajectory(:,3) + a)) * (2*pi) / C_;
            rot_samples(i,j,:,:) = reshape(trajectory,1,1,L_,[]);
        end
    end
    samples = rot_samples;
    %plot(reshape(samples(4,:,:,1),n_R_,L_)', reshape(samples(4,:,:,2),n_R_,L_)')
%ENDRM
end
