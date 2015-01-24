function turn = ccw(point_1, point_2, point_3)
%This function computes the determinant that gives the signed area of the
%triangle formed by point_1, point_2 and point_3. The function returns
%turn < 0 if the three points are in a counter-clockwise. Turn = 0 if they
%are collinear, and turn > 0 if they are clockwise.

x1=point_1(1,1); y1=point_1(2,1);
x2=point_2(1,1); y2=point_2(2,1);
x3=point_3(1,1); y3=point_3(2,1);

turn = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1);
