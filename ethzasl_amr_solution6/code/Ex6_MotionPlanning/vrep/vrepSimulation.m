%% V-REP Simulation Exercise 6: Motion Planning
% Tests the implemented Dijkstra's algorithm and the Dynamic Window Approach 
% within a V-Rep simulation framework.
% There is no need to change the code here, however you are free to play
% with the parameters.

% In order to run the simulation:
%   - Start V-Rep
%   - Load the scene file 'mooc_exercise.ttt'
%   - Start this Matlab script

more off;
startup;

%% Parameters setup

% if this is set to true, the path set constructed by the DWA is displayed
% during the simulation.
% ATTENTION: this makes the simulation much slower!
parameters.plotDwaPaths = false;

% this has to be set to the value chosen in V-REP
parameters.vrepTimeStep = 0.05;
 
% Define parameters for Dijkstra and Dynamic Window Approach
parameters.simTime = 5.0;
parameters.timestep = 2*parameters.vrepTimeStep; 
parameters.nVelSamples = 11;        % should be uneven
parameters.nOmegaSamples = 11;      % should be uneven
parameters.headingScoring = 0.1;
parameters.velocityScoring = 0.8;
parameters.obstacleDistanceScoring = 0.3;
parameters.maxVel = 0.2;
parameters.maxOmega = pi;
parameters.maxAcc = 1.0;
parameters.maxOmegaDot = pi;
parameters.plot = false;
parameters.connectivity = 8;
parameters.goalBrakingDistance = 0.5;
parameters.objectiveFcnSmoothingKernel = fspecial('gaussian', [3,3], 1.0);
parameters.globalPlanningOn = true; % use Global Dynamic Window Approach
parameters.safetyMargin = 0.1;      % make the robot bigger by this size to keep a safety margin due to delays and inaccuracies when modeling the dynamic constraints
parameters.vrepSteppedSimulation = true;

%% Initialize connection with V-Rep
connection = simulation_setup();

try

simulation_closeAllConnections(connection);  % make sure all connections are closed
connection = simulation_openConnection(connection, 0);

%% Get static data from V-Rep

simulation_start(connection);
bob_init(connection);
bodyDiameter = bob_getBodyDiameter(connection);
parameters.robotRadius = bodyDiameter/2.0 + parameters.safetyMargin;

wheelDiameter = bob_getWheelDiameter(connection);
wheelRadius = wheelDiameter/2.0;
interWheelDistance = bob_getInterWheelDistance(connection);

%STARTRM
scannerPoseWrtBob = bob_getScannerPose(connection);
%ENDRM

% Get the target pose
bob_setTargetGhostVisible(connection, true);
disp('Define the target pose in V-Rep and hit any key to continue...');
pause;
drawnow;
[goalPosition.x, goalPosition.y, dummy] = bob_getTargetGhostPose(connection);
tmp = goalPosition.x;                 % Transform the robot state to map (image) coordinates
goalPosition.x = -goalPosition.y;
goalPosition.y = tmp;
    
% Get bob pose
[robotState.x, robotState.y, robotState.heading] = bob_getPose(connection);

% Get global map
globalMapImg = bob_getMap(connection); % image has 512x512 values, corresponding to the 5x5 m^2 terrain
mapSize = size(globalMapImg);
mapOrigin = [-2.5, -2.5];       % TODO: get origin from V-Rep
mapResolution = 5.0/mapSize(1); % TODO: get from V-Rep

% Subsample map for faster Dijkstra to a resolution of 128x128
subsamplingFactor = 1/ceil(max(size(globalMapImg))/128);
globalMapImgSubsampled = imresize(globalMapImg, subsamplingFactor);
globalMapImgSubsampled(globalMapImgSubsampled > 0) = 1;
globalMapImgSubsampled(globalMapImgSubsampled <=0) = 0;
mapResolution = mapResolution/subsamplingFactor;
mapSize = size(globalMapImgSubsampled);
globalMap = createMap(mapOrigin, mapResolution, globalMapImgSubsampled);

%% Compute Dijkstra distance map
disp('Computing Dijkstra ...');
startIdx = worldToMap(globalMap.origin, globalMap.resolution, [robotState.x, robotState.y]);
goalIdx = worldToMap(globalMap.origin, globalMap.resolution, [goalPosition.x, goalPosition.y]);
[ costs, costGradientDirection ] = ...
    dijkstra( globalMap.data, goalIdx, parameters );
costGradientDirectionMap = createMap(globalMap.origin, globalMap.resolution, costGradientDirection);

%% DWA control in the loop
% Note: If you did not change the target pose, you do not have to run
% Dijkstra's algorithm again. You can simply execute this cell.

bob_setTargetGhostVisible(connection, true);
simulation_start(connection);

plotCosts(costs, [mapOrigin(2), mapOrigin(2)+mapResolution*mapSize(2)], ...
    [mapOrigin(1), mapOrigin(1)+mapResolution*mapSize(1)]);
drawnow;
hold on;

%STARTRM
% TODO: This should come from Exercise 2
%ENDRM
WheelSpeedsToVOmega = [wheelRadius/2 wheelRadius/2; ...
    -wheelRadius/interWheelDistance +wheelRadius/interWheelDistance];

[robotState.x, robotState.y, robotState.heading] = bob_getPose(connection);
goalDistance = hypot(robotState.x - goalPosition.x, robotState.y - goalPosition.y);

if parameters.vrepSteppedSimulation
    simulation_setStepped(connection, true);
end

nSpeedZeroCnt = 0; % Let's count the number of successive zero robot speeds to detect whether we are stuck
robotIsStuck = 0;
handles = [];

%STARTRM
% figure;
% hAxTest = axes;
% iter = 0;
%ENDRM

while goalDistance > 0.1

    % get robot pose (x,y,gamma), where positions are in meter, orientation in degrees (non-blocking function):
    [robotState.x robotState.y robotState.heading] = bob_getPose(connection);

    % Get the current v/omega of bob from V-REP
    [wheelSpeedLeftRadPerSec, wheelSpeedRightRadPerSec] = bob_getWheelSpeeds(connection);
    wheelSpeedsRadPerSec = [wheelSpeedLeftRadPerSec; wheelSpeedRightRadPerSec];
    
    % wheel speeds to v/omega
    vOmega = WheelSpeedsToVOmega*wheelSpeedsRadPerSec;
    robotState.vel = vOmega(1);
    robotState.omega = vOmega(2);

%STARTRM
    % TODO: do we want to create a map from the laser scans?

    % Get the latest laser scan. The scan is in the laser sensor frame.
    % The laser points will first be transformed into the global frame and
    % then entered into the local map.
%     [laserDataX laserDataY] = bob_getLaser(connection);
%     laserToGlobalHomogeneousTransformation2D = ...
%         [cos(robotState.heading + scannerPoseWrtBob(3)), -sin(robotState.heading + scannerPoseWrtBob(3)), robotState.x + scannerPoseWrtBob(1); ...
%          sin(robotState.heading + scannerPoseWrtBob(3)), cos(robotState.heading + scannerPoseWrtBob(3)), robotState.y + scannerPoseWrtBob(2); ...
%          0, 0, 1];
%     localMap = laserToMap(laserDataX, laserDataY, 0.1, laserToGlobalHomogeneousTransformation2D);
%ENDRM
    localMap = globalMap;

    % Transform the robot state to map (image) coordinates
    tmp = robotState.x;
    robotState.x = -robotState.y;
    robotState.y = tmp;
    robotState.heading = mod(robotState.heading + pi/2, 2*pi);

    % Call DWA
    [ v, omega, debug ] = dynamicWindowApproach( robotState, goalPosition, localMap, parameters, costGradientDirectionMap );
    
    % Let's enhance our planning a little bit: there is still a chance we
    % might get stuck. So let's detect that, and in the case of being
    % stuck, we start to rotate.
    if (abs(robotState.vel) < 1e-2 && abs(v) < 1e-2)
        nSpeedZeroCnt = nSpeedZeroCnt + 1;
    else
        nSpeedZeroCnt = 0;
    end
    if nSpeedZeroCnt > 20 % if N successive translational velocity commands are 0, we can assume we are stuck
        robotIsStuck = 1;
    else 
       robotIsStuck = 0; 
    end
    if robotIsStuck
        disp('Robot is stuck, starting to rotate');
        v = 0.0;
        omega = 0.5;
    end
   
     % v/omega to wheel speeds
    wheelSpeedsRadPerSec = WheelSpeedsToVOmega  \ [v; omega];
    % send command to V-REP
    bob_setWheelSpeeds(connection, wheelSpeedsRadPerSec(1), wheelSpeedsRadPerSec(2));
    % Did we reach the goal?
    goalDistance = hypot(robotState.x - goalPosition.x, robotState.y - goalPosition.y);

    % Debugging
    % plot the robot position
    if parameters.plotDwaPaths
        plot(robotState.y, robotState.x, 'og', 'MarkerFaceColor', 'g');
        % plot the trajectory set
        if debug.valid
            handles = plotTrajectories(debug, handles);
        end   
        drawnow;
    end

    % Trigger simulation step in V-rep
    if parameters.vrepSteppedSimulation
        nSimSteps = ceil(parameters.timestep/parameters.vrepTimeStep);
        nStepsDone = 0;
        while nStepsDone < nSimSteps 
            nStepsDone = nStepsDone + 1;
            simulation_triggerStep(connection); 
%STARTRM
%             iter = iter + 1;
%             [leftSpeed, rightSpeed] = bob_getWheelSpeeds(connection);
%             plot(hAxTest, iter, leftSpeed, 'b.', iter, wheelSpeedsRadPerSec(1), 'r.');
%             legend(hAxTest, 'true speed', 'commanded speed');
%             ylabel(hAxTest, 'wheel speed rad/s');
%             drawnow;
%             hold on;
%ENDRM
        end
    end
    
end

% Bring Bob to standstill
bob_setWheelSpeeds(connection, 0.0, 0.0);

if parameters.vrepSteppedSimulation
    simulation_setStepped(connection, false);
end

% Comment this in again if you want to stop the simulation automatically at
% the end
% simulation_stop(connection);
% simulation_closeConnection(connection);

catch exception
    simulation_closeConnection(connection);
    rethrow(exception);
end

disp('Simulation ended');
