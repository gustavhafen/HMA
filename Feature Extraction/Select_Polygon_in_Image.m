function polygon = Select_Polygon_in_Image(Image)
%function that returns a polygon drawed on image by the user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@parameters:
%image := the image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return:
%the polygon vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %handle to figure (to close it later)
    t = figure;
    imshow(uint8(Image)); 
    h = impoly; %let you create a polygon on the image
    wait(h); %to let select the points. Finish with double click
    polygon = getPosition(h); %the selected points as a matrix
    close(t);