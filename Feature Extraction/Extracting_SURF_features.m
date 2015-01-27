function Features = Extracting_SURF_features(img, Options)

img=im2double(img);
% Get the Key Points
Surf = OpenSurf(img,Options);
K = [ [Surf.x]; [Surf.y]; [Surf.scale]; [Surf.orientation] ];
D = [Surf.descriptor];

Features.Keypoints = double(K);
Features.Descriptors = double(D);
