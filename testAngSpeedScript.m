% script tests for angular speed function written for the
% plotHeadAngSpeed.m script. 
% Script tests for angular speed as calculated by trajectories representing three circles with increasing radii, each
% with three different numbers of coordinates. 
% Circles with different radii, so long as numbers of coordinates along the circumference are equal, should give identical angular speed. 
% The fewer coordinates there are along the circumference, the higher the
% angular speed should be.
% As a negative control, straight lines should give angular speed of zero. 
% script also tests for the effect of smoothing

% set parameters
smoothing = false;
frameRate = 9;


% generate circles
r=1; % radius
C=[1 1];
theta=0:2*pi/360:2*pi; % the angle
circle1=r*[cos(theta')+C(1) sin(theta')+C(2)]; % 360 points for the circle
circle2 = circle360(1:2:end,:); %180 points for the circle
circle3 = circle360(1:3:end,:); %120 points for the circle

% generate bigger circles
r=2; % radius
C=[1 1];
theta=0:2*pi/360:2*pi; % the angle
circle4=r*[cos(theta')+C(1) sin(theta')+C(2)]; % 360 points for the circle
circle5 = circle4(1:2:end,:); %180 points for the circle
circle6 = circle4(1:3:end,:); %120 points for the circle

% generate even bigger circles
r=3; % radius
C=[1 1];
theta=0:2*pi/360:2*pi; % the angle
circle7=r*[cos(theta')+C(1) sin(theta')+C(2)]; % 360 points for the circle
circle8 = circle7(1:2:end,:); %180 points for the circle
circle9 = circle7(1:3:end,:); %120 points for the circle


% generate a line
line360 = [1:360; 1:360]'; %360 points for a straight line
line180 = line360(1:2:end,:); %180 points for a straight line
line120 = line360(1:3:end,:); %120 points for a straight line

trajList = {circle1,circle2,circle3,circle4,circle5,circle6,circle7,circle8,circle9,line360,line180,line120};
angDiff = cell(length(trajList),1);
angSpeed = NaN(length(trajList),1);

% loop through each trajectory
for trajCtr = 1:length(trajList)
    traj = trajList{trajCtr};
    % get xy coordinates
    xcoords = traj(:,1)';
    ycoords = traj(:,2)';
    % calculate angles
    [angleArray,meanAngles] = makeAngleArray(xcoords,ycoords);
    angles= angleArray+meanAngles;
    % get angle difference over 1 second (over 1 second rather than between each frame to implement smoothing)
    if ~smoothing
        angleDiff = angles(2:end) - angles(1:end-1); % no smoothing
    else
        smoothFactor = frameRate+1;
        angleDiff = ((angles(smoothFactor:end) - angles(1:end-smoothFactor+1)))/frameRate;
    end
    % calculate total smoothed head angle change per second
    angSpeed(trajCtr) = abs(nansum(angleDiff)/length(angles)*frameRate);
end

% display values
angSpeed