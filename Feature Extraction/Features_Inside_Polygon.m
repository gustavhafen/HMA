function features = Features_Inside_Polygon(Features, Polygon)
%Only keypoints lying inside a polygon previously defined  by the user.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@paramters
%Features := An object containing the fields Keypoints and Descriptors.
%polygon     := polygon vertices (Matrix Nx2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@return:
%features := keypoints located inside the polygon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract keypoints

points = Features.Keypoints(1:2,:); %keypoints localizations
indexes = inpolygon(points(1,:),points(2,:),Polygon(:,1),Polygon(:,2)); %check if points fall inside the polygon
%% Packing
features.Keypoints = Features.Keypoints(:,indexes);
features.Descriptors = Features.Descriptors(:,indexes);