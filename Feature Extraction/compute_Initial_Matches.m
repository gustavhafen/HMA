function [polygon, Features, Initial_Matches, Images] = compute_Initial_Matches(image1, image2, params, polygon)
image1 = uint8(image1);
image2 = uint8(image2);
Images.Training = image1;
Images.Query = image2;
gray_Img1 = rgb2gray(image1);
gray_Img2 = rgb2gray(image2);
if nargin < 4,
    polygon = Select_Polygon_in_Image(image1);
end

Params_Features = Set_Param_Features(params.feature);

if strcmp(params.feature, 'SIFT'),
    Features1 = Extracting_SIFT_features(single(gray_Img1),Params_Features);
    Features2 = Extracting_SIFT_features(single(gray_Img2),Params_Features);
else if strcmp(params.feature, 'SURF'),
        Features1 = Extracting_SURF_features(gray_Img1, Params_Features.Options);
        Features2 = Extracting_SURF_features(gray_Img2, Params_Features.Options);
    else
        Features1 = Extracting_ASIFT_Features(gray_Img1, Params_Features);
        Features2 = Extracting_ASIFT_Features(gray_Img2, Params_Features);
    end
end


if ~isempty(polygon),
    Features1= Features_Inside_Polygon(Features1, polygon);
end

Initial_Matches = Set_Initial_Matches(Features1, Features2, params);
Features.Training = Features1;
Features.Query = Features2;
Images.crop_query = [];
Images.polygon_query_gt = [];
