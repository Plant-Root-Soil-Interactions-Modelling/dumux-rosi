RootRadius = 0.0002;
SoilRadius = 0.015;
Point(1) = {0,0,0};

Point(2) = {SoilRadius,0,0};
Point(3) = {0,SoilRadius,0};
Point(4) = {-SoilRadius,0,0};
Point(5) = {0,-SoilRadius,0};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Point(6) = {RootRadius,0,0};
Point(7) = {0,RootRadius,0};
Point(8) = {-RootRadius,0,0};
Point(9) = {0,-RootRadius,0};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};

Line(9) = {2,6};
Line(10) = {3,7};
Line(11) = {4,8};
Line(12) = {5,9};

Transfinite Line{5,6,7,8} = Ceil(3) Using Progression 1;
Transfinite Line{1,2,3,4} = Ceil(3) Using Progression 1;
Transfinite Line{-9,-10,-11,-12} = Ceil(30) Using Progression 1.2;

Line Loop(13) = {1,10,-5,-9};
Line Loop(14) = {2,11,-6,-10};
Line Loop(15) = {3,12,-7,-11};
Line Loop(16) = {4,9,-8,-12};
Line Loop(17) = {5,6,7,8};

Plane Surface(18) = {13};
Plane Surface(19) = {14};
Plane Surface(20) = {15};
Plane Surface(21) = {16};
Plane Surface(22) = {17};
Transfinite Surface {18,19,20,21,22};
Recombine Surface {18,19,20,21,22};

Extrude {0,0,-1} {
	  Surface{18,19,20,21,22}; Layers{1}; Recombine;
	}
Coherence;


Delete {
  Point{11};
}
Delete {
  Point{11};
}
Delete {
  Point{1};
}
Delete {
  Point{1};
}
Delete {
  Point{11};
}
Delete {
  Point{11, 1};
}
Delete {
  Point{1};
}
Delete {
  Point{11};
}
