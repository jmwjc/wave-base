
a = 205;
n = 40;

Point(1) = { -a,  -a, 0.0, a/n};
Point(2) = {  a,  -a, 0.0, a/n};
Point(3) = {  a,   a, 0.0, a/n};
Point(4) = { -a,   a, 0.0, a/n};
Point(5) = {  0,  -a, 0.0, a/n};

Line(1) = {1,5};
Line(2) = {5,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line(5) = {4,1};

Curve Loop(1) = {1,2,3,4,5};

Plane Surface(1) = {1};

Physical Curve("Γ") = {3,4,5};
Physical Surface("Ω") = {1};
Physical Point("Γᵗ") = {5};

// Transfinite Curve{1,2} = (n+1)/2;
// Transfinite Curve{3,4,5} = n;
// Transfinite Surface{1};
Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;