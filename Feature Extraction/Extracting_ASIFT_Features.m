function Features = Extracting_ASIFT_Features(img, Param)
%% Acquiring data %%%
    imwrite(img,'tmp.png');
    imwrite(zeros(100,100,3),'tmp2.png');    
    options = ' tmp.png tmp2.png tmp3 tmp4 tmp5 kpts.txt tmp6 [0]';
    unix([Param.exec_str  ' ' options]);
    loaded_keypoints = importdata('kpts.txt');
    num = loaded_keypoints(1,1);
    loaded_keypoints(1,:)=[];
    loaded_keypoints= reshape(loaded_keypoints', 132,num )';
    
    %Keypoints
    Features.Keypoints = loaded_keypoints(:,1:4)';
    %Descriptors
    Features.Descriptors = loaded_keypoints(:,5:end)';
