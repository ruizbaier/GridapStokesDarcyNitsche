/*************************************************************************\
 ComplexChannel.geo  - geometry script for a curved channel with 
                       three obstacles.
This is to be used with the Gmsh program
http://gmsh.info/
\*************************************************************************/

Mesh.ScalingFactor = 1.0;

// This value gives the global element size factor (lower -> finer mesh)
Mesh.CharacteristicLengthFactor = 0.02;

// This fixes the coloring to element type (better visibility)
Mesh.ColorCarousel = 0;

// Points are defined by {x,y,z,fac}, where fac gives a 
// relativ local refinement factor.

// control points for the channel
Point(1) = {-2, 1, 0, 1};
Point(2) = {-2, 0, 0, 1};
Point(3) = {0, 0, 0, 1};
Point(4) = {1.5, -.87, 0, 1};
Point(5) = {2.4, 0, 0, 1};
Point(6) = {0.1, 0.5, 0, 1};
Point(7) = {2, 1, 0, 1};
Point(8) = {0.5, 1, 0, 1};




//control points obstacle 1
Point(11) = {-1.4, 0.35, 0, 0.3};
Point(12) = {-1, 0.5, 0, 0.3};
Point(13) = {-0.6, 0.35, 0, 0.3};
Point(14) = {-0.6, 0.6, 0, 0.3};
Point(15) = {-1, 0.75, 0, 0.3};
Point(16) = {-1.4, 0.75, 0, 0.3};

t1 = 0.55; t2 =0.55; 
//control points obstacle 2
Point(21) = {-0.1+t1, -0.2+t2, 0, 0.3};
Point(22) = {0.2+t1, -0.2+t2, 0, 0.3};
Point(23) = {0.3+t1, -0.07+t2, 0, 0.3};
Point(24) = {0.1+t1, 0.05+t2, 0, 0.3};
Point(25) = {-0.1+t1, 0.2+t2, 0, 0.3};
Point(26) = {-0.3+t1, 0.1+t2, 0, 0.3};

v1 = 1.05; v2 = 0.6;

//control points obstacle 3
Point(31) = {0.55+v1, -0.8+v2, 0, 0.3};
Point(32) = {0.25+v1, -1.+v2, 0, 0.3};
Point(33) = {0.6+v1, -1.1+v2, 0, 0.3};
Point(34) = {0.8+v1, -0.9+v2, 0, 0.3};
Point(35) = {0.7+v1, -0.55+v2, 0, 0.3};
Point(36) = {0.5+v1, -0.65+v2, 0, 0.3};

// Lines are constructed from the points. Bsplines take control 
// points. Careful: lines follow a direction.

// channel shape
Line(1) = {1, 2};
BSpline(2) = {2, 3, 8, 4};
Line(3) = {4, 5};
BSpline(4) = {5, 6, 7, 1};

// Extra line
Point(9) = {-2, -0.2, 0, 1};
Point(10) = {1.37, -.99, 0, 1};
Point(83) = {-0., -0.3, 0, 1};
Point(88) = {0.4, 1, 0, 1};

Line(92) = {2, 9};
Line(93) = {4, 10};
BSpline(94) = {9, 83, 88, 10};



// obstacle 1
BSpline(5) = {11, 12, 13, 14, 15, 16, 11};

// obstacle 2
BSpline(6) = {21, 22, 23, 24, 25, 26, 21};

// obstacle 3
BSpline(7) = {31, 32, 33, 34, 35, 36, 31};

// Polygons and more general line loops are constructed from simple 
// lines. A minus sign substracts the object.

Line Loop(8) = {1, 2, 3, 4};
Line Loop(9) = {5};
Line Loop(10) = {6};
Line Loop(11) = {7};

// Surfaces are build similarly, second entries are considered holes.

Plane Surface(12) = {8, 9, 10, 11};

// Finally, we assign labels (10,20,30,40) to the lines or group 
// of lines. Gmsh seems to require a label for the surface, too.




Physical Point("inF", 10) = {1,2};
Physical Line("inF", 10) = {1};    // fluid inflow 

Physical Line("outF", 20) = {3};    // fluid outlow

Physical Point("wallF", 30) = {5}; 
Physical Line("wallF", 30) = {4};    // fluid wall

Physical Line("obsF", 40) = {5,6,7};    // fluid obstacles

Physical Point("interf",31) = {4}; 
Physical Line("interf",31) = {2}; // interface

//Physical Point("inP",41) = {9,2}; 
Physical Line("inP",41) = {92}; // porous inflow

Physical Line("wallP",42) = {94}; // porous wall

Physical Line("outP",43) = {93}; // porous outflow


Curve Loop(12) = {94, -93, -2, 92};
Plane Surface(44) = {12};


Physical Surface("fluid", 101) = {12};
Physical Surface("porous", 102) = {44};