
a = 205;
n = 5;

Point(1) = { -a,  -a, 0.0};
Point(2) = {  a,  -a, 0.0};
Point(3) = {  a,   a, 0.0};
Point(4) = { -a,   a, 0.0};
Point(5) = {  0,  -a, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Physical Curve("Γ") = {1,2,3,4};
Physical Surface("Ω") = {1};
Physical Point("Γ") = {5};

Transfinite Curve{1,2,3,4} = n+1;
Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;