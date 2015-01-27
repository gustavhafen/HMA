function Features = Extracting_SIFT_features(img, Params)

[K, D] = vl_sift(img, 'PeakThresh', Params.peak, 'EdgeThresh', Params.edge);

Features.Keypoints = double(K);
Features.Descriptors = double(D);