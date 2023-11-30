function centroid = computeCentreOfPressure(profile,cp)
    % Force on a panel
    Fi = [-cp.*profile.panels.normal(:,1).*profile.panels.length, ...
          -cp.*profile.panels.normal(:,2).*profile.panels.length, ...
           zeros(size(profile.panels.normal(:,1)))];
    % Horizontal force 
    Fx = Fi; Fx(:,2)=0; Fxtot = sum(Fx);
    % Vertical force
    Fy = Fi; Fy(:,1)=0; Fytot = sum(Fy);
    % arm of the i-th force around the origin
    ri = zeros(size(Fi)); 
    ri(:,1) = profile.panels.centre.X;
    ri(:,2) = profile.panels.centre.Y;
    ri(:,3) = 0;
    % horizontal arm
    rx = ri; rx(:,2) = 0;
    % vertical arm
    ry = ri; ry(:,2) = 0;

    % Moment of the vertical component of the i-th force
    Mx = zeros(size(Fi));
    for i = 1:numel(Mx(:,1))
        Mx(i,:) = cross(Fy(i,:),rx(i,:));
    end
    % Moment of the horizontal compoenent of the i-th force
    My = zeros(size(Fi));
    for i = 1:numel(My(:,1))
        My(i,:) = cross(Fx(i,:),ry(i,:));
    end
    Mxtot = sum(Mx); Mytot = sum(My);
    
    % Centre of pressure
    centroid.x = -Mxtot(3)/Fytot(2);
    centroid.y =  Mytot(3)/Fxtot(1);
end