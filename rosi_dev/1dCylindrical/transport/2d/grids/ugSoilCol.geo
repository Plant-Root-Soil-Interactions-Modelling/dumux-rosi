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

Point(10) = {4*RootRadius,0,0};
Point(11) = {0,4*RootRadius,0};
Point(12) = {-4*RootRadius,0,0};
Point(13) = {0,-4*RootRadius,0};
Circle(9) = {10,1,11};
Circle(10) = {11,1,12};
Circle(11) = {12,1,13};
Circle(12) = {13,1,10};


Transfinite Line{1,2,3,4,5,6,7,8,9,10,11,12} = Ceil(5) Using Progression 1;

Line Loop(1) = {5,6,7,8};
Line Loop(2) = {1,2,3,4};
Line Loop(3) = {9,10,11,12};

Plane Surface(1) = {2,3};
Plane Surface(2) = {1,3};


/*
Extrude {0,0,0.03} {
	  Surface{18,19,20,21,22}; Layers{1}; Recombine;
	}
Coherence;


Line Loop(9) = {1,2,3,4,5,6,7,8};
Plane Surface(10) = {9};
Transfinite Surface {10};

Extrude {0,0,0.03} {
  Surface{10};
  Layers{1};
  Recombine;
}

*/


Physical Surface(19) = {18};
