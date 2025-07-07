// mesh3d.geo  —  fully manual 3D extrusion of a 2D channel with 10 circular obstacles

// -----------------------------------------------------------------------------
// 0. Parameters & global mesh controls
// -----------------------------------------------------------------------------

// extrusion height
z0 = 0.8;

// small mesh size around obstacles
lc_obst = 0.5;

// Global mesh controls
Mesh.ScalingFactor              = 1.0;
Mesh.CharacteristicLengthFactor = 0.4;    // smaller → finer mesh globally
Mesh.ColorCarousel              = 0;

// -----------------------------------------------------------------------------
// 1. Original 2D geometry at z=0 (bottom layer)
// -----------------------------------------------------------------------------

// 1.1 Channel boundary points (z=0)
Point(1)  = { 0.0, 0.0, 0.0, 1.0 };
Point(2)  = { 0.0, 0.5, 0.0, 1.0 };
Point(3)  = {-0.2, 0.5, 0.0, 1.0 };
Point(4)  = {-0.2, 1.5, 0.0, 1.0 };
Point(5)  = { 0.0, 1.5, 0.0, 1.0 };
Point(6)  = { 0.0, 2.0, 0.0, 1.0 };
Point(7)  = { 3.2, 2.0, 0.0, 1.0 };
Point(8)  = { 3.2, 1.8, 0.0, 1.0 };
Point(9)  = { 4.2, 1.8, 0.0, 1.0 };
Point(10) = { 4.2, 1.3, 0.0, 1.0 };
Point(11) = { 3.2, 1.3, 0.0, 1.0 };
Point(12) = { 3.2, 0.7, 0.0, 1.0 };
Point(13) = { 4.2, 0.7, 0.0, 1.0 };
Point(14) = { 4.2, 0.2, 0.0, 1.0 };
Point(15) = { 3.2, 0.2, 0.0, 1.0 };
Point(16) = { 3.2, 0.0, 0.0, 1.0 };

// 1.2 Lines between channel points
Line(1)  = {  1,  2 };
Line(2)  = {  2,  3 };
Line(3)  = {  3,  4 };
Line(4)  = {  4,  5 };
Line(5)  = {  5,  6 };
Line(6)  = {  6,  7 };
Line(7)  = {  7,  8 };
Line(8)  = {  8,  9 };
Line(9)  = {  9, 10 };
Line(10) = { 10, 11 };
Line(11) = { 11, 12 };
Line(12) = { 12, 13 };
Line(13) = { 13, 14 };
Line(14) = { 14, 15 };
Line(15) = { 15, 16 };
Line(16) = { 16,  1 };

// 1.3 Outer loop at z=0
Curve Loop(1) = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };

// -----------------------------------------------------------------------------
// 2. Circular obstacles (exact circles via Arcs) at z=0
// -----------------------------------------------------------------------------

// Obstacle 1: center (0.7, 1.5), r = 0.2
Point(100) = { 0.7,      1.5,      0.0, lc_obst };
Point(101) = { 0.7 + 0.2, 1.5,      0.0, lc_obst };
Point(102) = { 0.7,      1.5 + 0.2, 0.0, lc_obst };
Point(103) = { 0.7 - 0.2, 1.5,      0.0, lc_obst };
Point(104) = { 0.7,      1.5 - 0.2, 0.0, lc_obst };
Circle(17) = { 101, 100, 102 };
Circle(18) = { 102, 100, 103 };
Circle(19) = { 103, 100, 104 };
Circle(20) = { 104, 100, 101 };
Curve Loop(2) = { 17, 18, 19, 20 };

// Obstacle 2: center (0.7, 1.0), r = 0.2
Point(105) = { 0.7,      1.0,      0.0, lc_obst };
Point(106) = { 0.7 + 0.2, 1.0,      0.0, lc_obst };
Point(107) = { 0.7,      1.0 + 0.2, 0.0, lc_obst };
Point(108) = { 0.7 - 0.2, 1.0,      0.0, lc_obst };
Point(109) = { 0.7,      1.0 - 0.2, 0.0, lc_obst };
Circle(21) = { 106, 105, 107 };
Circle(22) = { 107, 105, 108 };
Circle(23) = { 108, 105, 109 };
Circle(24) = { 109, 105, 106 };
Curve Loop(3) = { 21, 22, 23, 24 };

// Obstacle 3: center (0.7, 0.5), r = 0.2
Point(110) = { 0.7,      0.5,      0.0, lc_obst };
Point(111) = { 0.7 + 0.2, 0.5,      0.0, lc_obst };
Point(112) = { 0.7,      0.5 + 0.2, 0.0, lc_obst };
Point(113) = { 0.7 - 0.2, 0.5,      0.0, lc_obst };
Point(114) = { 0.7,      0.5 - 0.2, 0.0, lc_obst };
Circle(25) = { 111, 110, 112 };
Circle(26) = { 112, 110, 113 };
Circle(27) = { 113, 110, 114 };
Circle(28) = { 114, 110, 111 };
Curve Loop(4) = { 25, 26, 27, 28 };

// Obstacle 4: center (2.1, 1.5), r = 0.2
Point(115) = { 2.1,      1.5,      0.0, lc_obst };
Point(116) = { 2.1 + 0.2, 1.5,      0.0, lc_obst };
Point(117) = { 2.1,      1.5 + 0.2, 0.0, lc_obst };
Point(118) = { 2.1 - 0.2, 1.5,      0.0, lc_obst };
Point(119) = { 2.1,      1.5 - 0.2, 0.0, lc_obst };
Circle(29) = { 116, 115, 117 };
Circle(30) = { 117, 115, 118 };
Circle(31) = { 118, 115, 119 };
Circle(32) = { 119, 115, 116 };
Curve Loop(5) = { 29, 30, 31, 32 };

// Obstacle 5: center (2.1, 1.0), r = 0.2
Point(120) = { 2.1,      1.0,      0.0, lc_obst };
Point(121) = { 2.1 + 0.2, 1.0,      0.0, lc_obst };
Point(122) = { 2.1,      1.0 + 0.2, 0.0, lc_obst };
Point(123) = { 2.1 - 0.2, 1.0,      0.0, lc_obst };
Point(124) = { 2.1,      1.0 - 0.2, 0.0, lc_obst };
Circle(33) = { 121, 120, 122 };
Circle(34) = { 122, 120, 123 };
Circle(35) = { 123, 120, 124 };
Circle(36) = { 124, 120, 121 };
Curve Loop(6) = { 33, 34, 35, 36 };

// Obstacle 6: center (2.1, 0.5), r = 0.2
Point(125) = { 2.1,      0.5,      0.0, lc_obst };
Point(126) = { 2.1 + 0.2, 0.5,      0.0, lc_obst };
Point(127) = { 2.1,      0.5 + 0.2, 0.0, lc_obst };
Point(128) = { 2.1 - 0.2, 0.5,      0.0, lc_obst };
Point(129) = { 2.1,      0.5 - 0.2, 0.0, lc_obst };
Circle(37) = { 126, 125, 127 };
Circle(38) = { 127, 125, 128 };
Circle(39) = { 128, 125, 129 };
Circle(40) = { 129, 125, 126 };
Curve Loop(7) = { 37, 38, 39, 40 };

// Obstacle 7: center (1.4, 1.645), r = 0.15
Point(130) = { 1.4,      1.645,     0.0, lc_obst };
Point(131) = { 1.4 + 0.15, 1.645,     0.0, lc_obst };
Point(132) = { 1.4,      1.645 + 0.15, 0.0, lc_obst };
Point(133) = { 1.4 - 0.15, 1.645,     0.0, lc_obst };
Point(134) = { 1.4,      1.645 - 0.15, 0.0, lc_obst };
Circle(41) = { 131, 130, 132 };
Circle(42) = { 132, 130, 133 };
Circle(43) = { 133, 130, 134 };
Circle(44) = { 134, 130, 131 };
Curve Loop(8) = { 41, 42, 43, 44 };

// Obstacle 8: center (1.4, 1.215), r = 0.15
Point(135) = { 1.4,      1.215,     0.0, lc_obst };
Point(136) = { 1.4 + 0.15, 1.215,     0.0, lc_obst };
Point(137) = { 1.4,      1.215 + 0.15, 0.0, lc_obst };
Point(138) = { 1.4 - 0.15, 1.215,     0.0, lc_obst };
Point(139) = { 1.4,      1.215 - 0.15, 0.0, lc_obst };
Circle(45) = { 136, 135, 137 };
Circle(46) = { 137, 135, 138 };
Circle(47) = { 138, 135, 139 };
Circle(48) = { 139, 135, 136 };
Curve Loop(9) = { 45, 46, 47, 48 };

// Obstacle 9: center (1.4, 0.785), r = 0.15
Point(140) = { 1.4,      0.785,     0.0, lc_obst };
Point(141) = { 1.4 + 0.15, 0.785,     0.0, lc_obst };
Point(142) = { 1.4,      0.785 + 0.15, 0.0, lc_obst };
Point(143) = { 1.4 - 0.15, 0.785,     0.0, lc_obst };
Point(144) = { 1.4,      0.785 - 0.15, 0.0, lc_obst };
Circle(49) = { 141, 140, 142 };
Circle(50) = { 142, 140, 143 };
Circle(51) = { 143, 140, 144 };
Circle(52) = { 144, 140, 141 };
Curve Loop(10) = { 49, 50, 51, 52 };

// Obstacle 10: center (1.4, 0.355), r = 0.15
Point(145) = { 1.4,      0.355,     0.0, lc_obst };
Point(146) = { 1.4 + 0.15, 0.355,     0.0, lc_obst };
Point(147) = { 1.4,      0.355 + 0.15, 0.0, lc_obst };
Point(148) = { 1.4 - 0.15, 0.355,     0.0, lc_obst };
Point(149) = { 1.4,      0.355 - 0.15, 0.0, lc_obst };
Circle(53) = { 146, 145, 147 };
Circle(54) = { 147, 145, 148 };
Circle(55) = { 148, 145, 149 };
Circle(56) = { 149, 145, 146 };
Curve Loop(11) = { 53, 54, 55, 56 };

// 1.4 **Plane Surfaces** (bottom):
//     - Surface(1) : channel floor minus ten circular holes
Plane Surface(1)  = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

// - Surfaces (2..11): the ten obstacle bottoms (each a standalone disk at z=0)
Plane Surface(2)  = { 2 };
Plane Surface(3)  = { 3 };
Plane Surface(4)  = { 4 };
Plane Surface(5)  = { 5 };
Plane Surface(6)  = { 6 };
Plane Surface(7)  = { 7 };
Plane Surface(8)  = { 8 };
Plane Surface(9)  = { 9 };
Plane Surface(10) = { 10 };
Plane Surface(11) = { 11 };

// -----------------------------------------------------------------------------
// 2. Duplicate geometry at z = z0 (top layer)
// -----------------------------------------------------------------------------

// 2.1 Outer loop points at z = z0
Point(1001) = { 0.0, 0.0, z0, 1.0 };
Point(1002) = { 0.0, 0.5, z0, 1.0 };
Point(1003) = {-0.2, 0.5, z0, 1.0 };
Point(1004) = {-0.2, 1.5, z0, 1.0 };
Point(1005) = { 0.0, 1.5, z0, 1.0 };
Point(1006) = { 0.0, 2.0, z0, 1.0 };
Point(1007) = { 3.2, 2.0, z0, 1.0 };
Point(1008) = { 3.2, 1.8, z0, 1.0 };
Point(1009) = { 4.2, 1.8, z0, 1.0 };
Point(1010) = { 4.2, 1.3, z0, 1.0 };
Point(1011) = { 3.2, 1.3, z0, 1.0 };
Point(1012) = { 3.2, 0.7, z0, 1.0 };
Point(1013) = { 4.2, 0.7, z0, 1.0 };
Point(1014) = { 4.2, 0.2, z0, 1.0 };
Point(1015) = { 3.2, 0.2, z0, 1.0 };
Point(1016) = { 3.2, 0.0, z0, 1.0 };

// 2.2 Obstacle circle centers at z = z0
Point(1100) = { 0.7,    1.5,      z0, lc_obst };
Point(1101) = { 0.9,    1.5,      z0, lc_obst };
Point(1102) = { 0.7,    1.7,      z0, lc_obst };
Point(1103) = { 0.5,    1.5,      z0, lc_obst };
Point(1104) = { 0.7,    1.3,      z0, lc_obst };

Point(1105) = { 0.7,    1.0,      z0, lc_obst };
Point(1106) = { 0.9,    1.0,      z0, lc_obst };
Point(1107) = { 0.7,    1.2,      z0, lc_obst };
Point(1108) = { 0.5,    1.0,      z0, lc_obst };
Point(1109) = { 0.7,    0.8,      z0, lc_obst };

Point(1110) = { 0.7,    0.5,      z0, lc_obst };
Point(1111) = { 0.9,    0.5,      z0, lc_obst };
Point(1112) = { 0.7,    0.7,      z0, lc_obst };
Point(1113) = { 0.5,    0.5,      z0, lc_obst };
Point(1114) = { 0.7,    0.3,      z0, lc_obst };

Point(1115) = { 2.1,    1.5,      z0, lc_obst };
Point(1116) = { 2.3,    1.5,      z0, lc_obst };
Point(1117) = { 2.1,    1.7,      z0, lc_obst };
Point(1118) = { 1.9,    1.5,      z0, lc_obst };
Point(1119) = { 2.1,    1.3,      z0, lc_obst };

Point(1120) = { 2.1,    1.0,      z0, lc_obst };
Point(1121) = { 2.3,    1.0,      z0, lc_obst };
Point(1122) = { 2.1,    1.2,      z0, lc_obst };
Point(1123) = { 1.9,    1.0,      z0, lc_obst };
Point(1124) = { 2.1,    0.8,      z0, lc_obst };

Point(1125) = { 2.1,    0.5,      z0, lc_obst };
Point(1126) = { 2.3,    0.5,      z0, lc_obst };
Point(1127) = { 2.1,    0.7,      z0, lc_obst };
Point(1128) = { 1.9,    0.5,      z0, lc_obst };
Point(1129) = { 2.1,    0.3,      z0, lc_obst };

Point(1130) = { 1.4,    1.645,    z0, lc_obst };
Point(1131) = { 1.55,   1.645,    z0, lc_obst };
Point(1132) = { 1.4,    1.795,    z0, lc_obst };
Point(1133) = { 1.25,   1.645,    z0, lc_obst };
Point(1134) = { 1.4,    1.495,    z0, lc_obst };

Point(1135) = { 1.4,    1.215,    z0, lc_obst };
Point(1136) = { 1.55,   1.215,    z0, lc_obst };
Point(1137) = { 1.4,    1.365,    z0, lc_obst };
Point(1138) = { 1.25,   1.215,    z0, lc_obst };
Point(1139) = { 1.4,    1.065,    z0, lc_obst };

Point(1140) = { 1.4,    0.785,    z0, lc_obst };
Point(1141) = { 1.55,   0.785,    z0, lc_obst };
Point(1142) = { 1.4,    0.935,    z0, lc_obst };
Point(1143) = { 1.25,   0.785,    z0, lc_obst };
Point(1144) = { 1.4,    0.635,    z0, lc_obst };

Point(1145) = { 1.4,    0.355,    z0, lc_obst };
Point(1146) = { 1.55,   0.355,    z0, lc_obst };
Point(1147) = { 1.4,    0.505,    z0, lc_obst };
Point(1148) = { 1.25,   0.355,    z0, lc_obst };
Point(1149) = { 1.4,    0.205,    z0, lc_obst };

// 2.3 Duplicate outer lines (connect bottom→top)
Line(2001) = { 1001,1002 };
Line(2002) = { 1002,1003 };
Line(2003) = { 1003,1004 };
Line(2004) = { 1004,1005 };
Line(2005) = { 1005,1006 };
Line(2006) = { 1006,1007 };
Line(2007) = { 1007,1008 };
Line(2008) = { 1008,1009 };
Line(2009) = { 1009,1010 };
Line(2010) = { 1010,1011 };
Line(2011) = { 1011,1012 };
Line(2012) = { 1012,1013 };
Line(2013) = { 1013,1014 };
Line(2014) = { 1014,1015 };
Line(2015) = { 1015,1016 };
Line(2016) = { 1016,1001 };

// 2.4 Duplicate circle arcs (connect bottom→top)
Circle(2017) = { 1101, 1100, 1102 };  // Obstacle 1 top arcs
Circle(2018) = { 1102, 1100, 1103 };
Circle(2019) = { 1103, 1100, 1104 };
Circle(2020) = { 1104, 1100, 1101 };

Circle(2021) = { 1106, 1105, 1107 };  // Obstacle 2 top arcs
Circle(2022) = { 1107, 1105, 1108 };
Circle(2023) = { 1108, 1105, 1109 };
Circle(2024) = { 1109, 1105, 1106 };

Circle(2025) = { 1111, 1110, 1112 };  // Obstacle 3 top arcs
Circle(2026) = { 1112, 1110, 1113 };
Circle(2027) = { 1113, 1110, 1114 };
Circle(2028) = { 1114, 1110, 1111 };

Circle(2029) = { 1116, 1115, 1117 };  // Obstacle 4 top arcs
Circle(2030) = { 1117, 1115, 1118 };
Circle(2031) = { 1118, 1115, 1119 };
Circle(2032) = { 1119, 1115, 1116 };

Circle(2033) = { 1121, 1120, 1122 };  // Obstacle 5 top arcs
Circle(2034) = { 1122, 1120, 1123 };
Circle(2035) = { 1123, 1120, 1124 };
Circle(2036) = { 1124, 1120, 1121 };

Circle(2037) = { 1126, 1125, 1127 };  // Obstacle 6 top arcs
Circle(2038) = { 1127, 1125, 1128 };
Circle(2039) = { 1128, 1125, 1129 };
Circle(2040) = { 1129, 1125, 1126 };

Circle(2041) = { 1131, 1130, 1132 };  // Obstacle 7 top arcs
Circle(2042) = { 1132, 1130, 1133 };
Circle(2043) = { 1133, 1130, 1134 };
Circle(2044) = { 1134, 1130, 1131 };

Circle(2045) = { 1136, 1135, 1137 };  // Obstacle 8 top arcs
Circle(2046) = { 1137, 1135, 1138 };
Circle(2047) = { 1138, 1135, 1139 };
Circle(2048) = { 1139, 1135, 1136 };

Circle(2049) = { 1141, 1140, 1142 };  // Obstacle 9 top arcs
Circle(2050) = { 1142, 1140, 1143 };
Circle(2051) = { 1143, 1140, 1144 };
Circle(2052) = { 1144, 1140, 1141 };

Circle(2053) = { 1146, 1145, 1147 };  // Obstacle 10 top arcs
Circle(2054) = { 1147, 1145, 1148 };
Circle(2055) = { 1148, 1145, 1149 };
Circle(2056) = { 1149, 1145, 1146 };

// 2.5 **cylinder top curve‐loops**  (z = z0)  
//      used by Plane Surface(102..111) for each cylinder
Curve Loop(102) = { 2017, 2018, 2019, 2020 };   // Obstacle 1 top
Curve Loop(103) = { 2021, 2022, 2023, 2024 };   // Obstacle 2 top
Curve Loop(104) = { 2025, 2026, 2027, 2028 };   // Obstacle 3 top
Curve Loop(105) = { 2029, 2030, 2031, 2032 };   // Obstacle 4 top
Curve Loop(106) = { 2033, 2034, 2035, 2036 };   // Obstacle 5 top
Curve Loop(107) = { 2037, 2038, 2039, 2040 };   // Obstacle 6 top
Curve Loop(108) = { 2041, 2042, 2043, 2044 };   // Obstacle 7 top
Curve Loop(109) = { 2045, 2046, 2047, 2048 };   // Obstacle 8 top
Curve Loop(110) = { 2049, 2050, 2051, 2052 };   // Obstacle 9 top
Curve Loop(111) = { 2053, 2054, 2055, 2056 };   // Obstacle 10 top

// 2.6 **cylinder top surfaces**  (z = z0)
Plane Surface(102) = { 102 };   // Obstacle 1 top
Plane Surface(103) = { 103 };   // Obstacle 2 top
Plane Surface(104) = { 104 };   // Obstacle 3 top
Plane Surface(105) = { 105 };   // Obstacle 4 top
Plane Surface(106) = { 106 };   // Obstacle 5 top
Plane Surface(107) = { 107 };   // Obstacle 6 top
Plane Surface(108) = { 108 };   // Obstacle 7 top
Plane Surface(109) = { 109 };   // Obstacle 8 top
Plane Surface(110) = { 110 };   // Obstacle 9 top
Plane Surface(111) = { 111 };   // Obstacle 10 top

// 2.7 **Fluid‐ceiling “hole” curve‐loops**  (z = z0)  
//      exactly duplicate the arcs, but assigned new loop IDs
Curve Loop(112) = { 2017, 2018, 2019, 2020 };   // hole for cylinder 1
Curve Loop(113) = { 2021, 2022, 2023, 2024 };   // hole for cylinder 2
Curve Loop(114) = { 2025, 2026, 2027, 2028 };   // hole for cylinder 3
Curve Loop(115) = { 2029, 2030, 2031, 2032 };   // hole for cylinder 4
Curve Loop(116) = { 2033, 2034, 2035, 2036 };   // hole for cylinder 5
Curve Loop(117) = { 2037, 2038, 2039, 2040 };   // hole for cylinder 6
Curve Loop(118) = { 2041, 2042, 2043, 2044 };   // hole for cylinder 7
Curve Loop(119) = { 2045, 2046, 2047, 2048 };   // hole for cylinder 8
Curve Loop(120) = { 2049, 2050, 2051, 2052 };   // hole for cylinder 9
Curve Loop(121) = { 2053, 2054, 2055, 2056 };   // hole for cylinder 10

// 2.8 **Fluid‐ceiling surface**  (z = z0)
//      one single surface that uses outer loop 101 minus holes 112..121
Curve Loop(101) = { 2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016 };
Plane Surface(101) = { 101, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121 };

// -----------------------------------------------------------------------------
// 3. Connector lines between bottom and top layers
// -----------------------------------------------------------------------------

// Connect each bottom point i to its duplicate 100i
Line(3001) = { 1,   1001 }; Line(3002) = { 2,   1002 };
Line(3003) = { 3,   1003 }; Line(3004) = { 4,   1004 };
Line(3005) = { 5,   1005 }; Line(3006) = { 6,   1006 };
Line(3007) = { 7,   1007 }; Line(3008) = { 8,   1008 };
Line(3009) = { 9,   1009 }; Line(3010) = {10,   1010 };
Line(3011) = {11,   1011 }; Line(3012) = {12,   1012 };
Line(3013) = {13,   1013 }; Line(3014) = {14,   1014 };
Line(3015) = {15,   1015 }; Line(3016) = {16,   1016 };

// Connect each obstacle bottom‐center to its top‐center
Line(3100) = { 100, 1100 }; Line(3101) = { 101, 1101 };
Line(3102) = { 102, 1102 }; Line(3103) = { 103, 1103 };
Line(3104) = { 104, 1104 }; Line(3105) = { 105, 1105 };
Line(3106) = { 106, 1106 }; Line(3107) = { 107, 1107 };
Line(3108) = { 108, 1108 }; Line(3109) = { 109, 1109 };
Line(3110) = { 110, 1110 }; Line(3111) = { 111, 1111 };
Line(3112) = { 112, 1112 }; Line(3113) = { 113, 1113 };
Line(3114) = { 114, 1114 }; Line(3115) = { 115, 1115 };
Line(3116) = { 116, 1116 }; Line(3117) = { 117, 1117 };
Line(3118) = { 118, 1118 }; Line(3119) = { 119, 1119 };
Line(3120) = { 120, 1120 }; Line(3121) = { 121, 1121 };
Line(3122) = { 122, 1122 }; Line(3123) = { 123, 1123 };
Line(3124) = { 124, 1124 }; Line(3125) = { 125, 1125 };
Line(3126) = { 126, 1126 }; Line(3127) = { 127, 1127 };
Line(3128) = { 128, 1128 }; Line(3129) = { 129, 1129 };
Line(3130) = { 130, 1130 }; Line(3131) = { 131, 1131 };
Line(3132) = { 132, 1132 }; Line(3133) = { 133, 1133 };
Line(3134) = { 134, 1134 }; Line(3135) = { 135, 1135 };
Line(3136) = { 136, 1136 }; Line(3137) = { 137, 1137 };
Line(3138) = { 138, 1138 }; Line(3139) = { 139, 1139 };
Line(3140) = { 140, 1140 }; Line(3141) = { 141, 1141 };
Line(3142) = { 142, 1142 }; Line(3143) = { 143, 1143 };
Line(3144) = { 144, 1144 }; Line(3145) = { 145, 1145 };
Line(3146) = { 146, 1146 }; Line(3147) = { 147, 1147 };
Line(3148) = { 148, 1148 }; Line(3149) = { 149, 1149 };

// -----------------------------------------------------------------------------
// 4. Side‐wall surfaces between bottom & top curves
// -----------------------------------------------------------------------------

// 4.1 Outer channel side‐walls (16 faces)
Curve Loop(4001) = {   1, 3002, -2001, -3001 };  Surface(4001) = {4001};
Curve Loop(4002) = {   2, 3003, -2002, -3002 };  Surface(4002) = {4002};
Curve Loop(4003) = {   3, 3004, -2003, -3003 };  Surface(4003) = {4003};
Curve Loop(4004) = {   4, 3005, -2004, -3004 };  Surface(4004) = {4004};
Curve Loop(4005) = {   5, 3006, -2005, -3005 };  Surface(4005) = {4005};
Curve Loop(4006) = {   6, 3007, -2006, -3006 };  Surface(4006) = {4006};
Curve Loop(4007) = {   7, 3008, -2007, -3007 };  Surface(4007) = {4007};
Curve Loop(4008) = {   8, 3009, -2008, -3008 };  Surface(4008) = {4008};
Curve Loop(4009) = {   9, 3010, -2009, -3009 };  Surface(4009) = {4009};
Curve Loop(4010) = {  10, 3011, -2010, -3010 };  Surface(4010) = {4010};
Curve Loop(4011) = {  11, 3012, -2011, -3011 };  Surface(4011) = {4011};
Curve Loop(4012) = {  12, 3013, -2012, -3012 };  Surface(4012) = {4012};
Curve Loop(4013) = {  13, 3014, -2013, -3013 };  Surface(4013) = {4013};
Curve Loop(4014) = {  14, 3015, -2014, -3014 };  Surface(4014) = {4014};
Curve Loop(4015) = {  15, 3016, -2015, -3015 };  Surface(4015) = {4015};
Curve Loop(4016) = {  16, 3001, -2016, -3016 };  Surface(4016) = {4016};

// 4.2 Obstacle 1 side‐walls (4 faces)
Curve Loop(4017) = {  17, 3102, -2017, -3101 };  Surface(4017) = {4017};
Curve Loop(4018) = {  18, 3103, -2018, -3102 };  Surface(4018) = {4018};
Curve Loop(4019) = {  19, 3104, -2019, -3103 };  Surface(4019) = {4019};
Curve Loop(4020) = {  20, 3101, -2020, -3104 };  Surface(4020) = {4020};

// 4.3 Obstacle 2 side‐walls
Curve Loop(4021) = {  21, 3107, -2021, -3106 };  Surface(4021) = {4021};
Curve Loop(4022) = {  22, 3108, -2022, -3107 };  Surface(4022) = {4022};
Curve Loop(4023) = {  23, 3109, -2023, -3108 };  Surface(4023) = {4023};
Curve Loop(4024) = {  24, 3106, -2024, -3109 };  Surface(4024) = {4024};

// 4.4 Obstacle 3 side‐walls
Curve Loop(4025) = {  25, 3112, -2025, -3111 };  Surface(4025) = {4025};
Curve Loop(4026) = {  26, 3113, -2026, -3112 };  Surface(4026) = {4026};
Curve Loop(4027) = {  27, 3114, -2027, -3113 };  Surface(4027) = {4027};
Curve Loop(4028) = {  28, 3111, -2028, -3114 };  Surface(4028) = {4028};

// 4.5 Obstacle 4 side‐walls
Curve Loop(4029) = {  29, 3117, -2029, -3116 };  Surface(4029) = {4029};
Curve Loop(4030) = {  30, 3118, -2030, -3117 };  Surface(4030) = {4030};
Curve Loop(4031) = {  31, 3119, -2031, -3118 };  Surface(4031) = {4031};
Curve Loop(4032) = {  32, 3116, -2032, -3119 };  Surface(4032) = {4032};

// 4.6 Obstacle 5 side‐walls
Curve Loop(4033) = {  33, 3122, -2033, -3121 };  Surface(4033) = {4033};
Curve Loop(4034) = {  34, 3123, -2034, -3122 };  Surface(4034) = {4034};
Curve Loop(4035) = {  35, 3124, -2035, -3123 };  Surface(4035) = {4035};
Curve Loop(4036) = {  36, 3121, -2036, -3124 };  Surface(4036) = {4036};

// 4.7 Obstacle 6 side‐walls
Curve Loop(4037) = {  37, 3127, -2037, -3126 };  Surface(4037) = {4037};
Curve Loop(4038) = {  38, 3128, -2038, -3127 };  Surface(4038) = {4038};
Curve Loop(4039) = {  39, 3129, -2039, -3128 };  Surface(4039) = {4039};
Curve Loop(4040) = {  40, 3126, -2040, -3129 };  Surface(4040) = {4040};

// 4.8 Obstacle 7 side‐walls
Curve Loop(4041) = {  41, 3132, -2041, -3131 };  Surface(4041) = {4041};
Curve Loop(4042) = {  42, 3133, -2042, -3132 };  Surface(4042) = {4042};
Curve Loop(4043) = {  43, 3134, -2043, -3133 };  Surface(4043) = {4043};
Curve Loop(4044) = {  44, 3131, -2044, -3134 };  Surface(4044) = {4044};

// 4.9 Obstacle 8 side‐walls
Curve Loop(4045) = {  45, 3137, -2045, -3136 };  Surface(4045) = {4045};
Curve Loop(4046) = {  46, 3138, -2046, -3137 };  Surface(4046) = {4046};
Curve Loop(4047) = {  47, 3139, -2047, -3138 };  Surface(4047) = {4047};
Curve Loop(4048) = {  48, 3136, -2048, -3139 };  Surface(4048) = {4048};

// 4.10 Obstacle 9 side‐walls
Curve Loop(4049) = {  49, 3142, -2049, -3141 };  Surface(4049) = {4049};
Curve Loop(4050) = {  50, 3143, -2050, -3142 };  Surface(4050) = {4050};
Curve Loop(4051) = {  51, 3144, -2051, -3143 };  Surface(4051) = {4051};
Curve Loop(4052) = {  52, 3141, -2052, -3144 };  Surface(4052) = {4052};

// 4.11 Obstacle 10 side‐walls
Curve Loop(4053) = {  53, 3147, -2053, -3146 };  Surface(4053) = {4053};
Curve Loop(4054) = {  54, 3148, -2054, -3147 };  Surface(4054) = {4054};
Curve Loop(4055) = {  55, 3149, -2055, -3148 };  Surface(4055) = {4055};
Curve Loop(4056) = {  56, 3146, -2056, -3149 };  Surface(4056) = {4056};

// -----------------------------------------------------------------------------
// 5. Define volumes for each of the ten cylinders (porous cores)
// -----------------------------------------------------------------------------

// Obstacle #1 
Surface Loop(5001) = {   2,  102, 4017, 4018, 4019, 4020 };
Volume(5001)     = { 5001 };

// Obstacle #2 
Surface Loop(5002) = {   3,  103, 4021, 4022, 4023, 4024 };
Volume(5002)     = { 5002 };

// Obstacle #3 
Surface Loop(5003) = {   4,  104, 4025, 4026, 4027, 4028 };
Volume(5003)     = { 5003 };

// Obstacle #4 
Surface Loop(5004) = {   5,  105, 4029, 4030, 4031, 4032 };
Volume(5004)     = { 5004 };

// Obstacle #5 
Surface Loop(5005) = {   6,  106, 4033, 4034, 4035, 4036 };
Volume(5005)     = { 5005 };

// Obstacle #6 
Surface Loop(5006) = {   7,  107, 4037, 4038, 4039, 4040 };
Volume(5006)     = { 5006 };

// Obstacle #7 
Surface Loop(5007) = {   8,  108, 4041, 4042, 4043, 4044 };
Volume(5007)     = { 5007 };

// Obstacle #8 
Surface Loop(5008) = {   9,  109, 4045, 4046, 4047, 4048 };
Volume(5008)     = { 5008 };

// Obstacle #9 
Surface Loop(5009) = {  10,  110, 4049, 4050, 4051, 4052 };
Volume(5009)     = { 5009 };

// Obstacle #10 
Surface Loop(5010) = {  11,  111, 4053, 4054, 4055, 4056 };
Volume(5010)     = { 5010 };

// -----------------------------------------------------------------------------
// 6. Assemble the Fluid volume (everything outside the obstacles, inside channel)
// -----------------------------------------------------------------------------

//    – bottom floor   = Surface(1)   (outer boundary minus holes 2..11)
//    – top ceiling    = Surface(101) (outer boundary minus holes 112..121)
//    – outer side walls = Surfaces(4001..4016) (positive orientation)
//    – obstacle side walls (4017..4056) appear with negative sign (shared with cylinder volumes)
Surface Loop(5100) = {
       1,        // bottom floor with holes
      101,       // top ceiling with holes
   4001,4002,4003,4004,4005,4006,4007,4008,
   4009,4010,4011,4012,4013,4014,4015,4016,
  -4017,-4018,-4019,-4020,    // cylinder 1 sides (negative)
  -4021,-4022,-4023,-4024,    // cylinder 2 sides
  -4025,-4026,-4027,-4028,    // cylinder 3 sides
  -4029,-4030,-4031,-4032,    // cylinder 4 sides
  -4033,-4034,-4035,-4036,    // cylinder 5 sides
  -4037,-4038,-4039,-4040,    // cylinder 6 sides
  -4041,-4042,-4043,-4044,    // cylinder 7 sides
  -4045,-4046,-4047,-4048,    // cylinder 8 sides
  -4049,-4050,-4051,-4052,    // cylinder 9 sides
  -4053,-4054,-4055,-4056     // cylinder 10 sides
};
Volume(5100) = { 5100 };

// -----------------------------------------------------------------------------
// 7. Tag individual surfaces as Physical Surfaces
//
// - "Inlet" : the upstream face (x=0 side)
// - "Outlet" : the downstream face (x=4.2 side of the right branch)
// - "LateralWalls" : all other channel side‐walls
// - "Interfaces" : all 10 cylindrical side‐walls
// - "CylinderTops" : the ten top disks (z=z₀)
// - "CylinderBottoms" : the ten bottom disks (z=0)
// - "ChannelTop" : the ceiling minus holes (surface 101)
// - "ChannelBottom" : the floor minus holes (surface 1)
// -----------------------------------------------------------------------------

// 7.1 Inlet and Outlet faces
Physical Surface("Inlet") = {  4003 }; // side‐face extruded from Line(3)
Physical Surface("Outlet") = { 4009, 4013 }; // side‐face extruded from Line(13)

// 7.2 All other channel side‐walls (lateral surfaces, excluding inlet/outlet)
Physical Surface("LateralWalls") = { 4001, 4002, 4004, 4005, 4006, 4007, 4008, 4010, 4011, 4012, 4014, 4015, 4016 };

// 7.3 All 10 cylindrical side‐walls (fluid/obstacle interfaces)
Physical Surface("Interf") = {
    4017, 4018, 4019, 4020,   // obstacle 1 side‐walls
    4021, 4022, 4023, 4024,   // obstacle 2 side‐walls
    4025, 4026, 4027, 4028,   // obstacle 3 side‐walls
    4029, 4030, 4031, 4032,   // obstacle 4 side‐walls
    4033, 4034, 4035, 4036,   // obstacle 5 side‐walls
    4037, 4038, 4039, 4040,   // obstacle 6 side‐walls
    4041, 4042, 4043, 4044,   // obstacle 7 side‐walls
    4045, 4046, 4047, 4048,   // obstacle 8 side‐walls
    4049, 4050, 4051, 4052,   // obstacle 9 side‐walls
    4053, 4054, 4055, 4056    // obstacle 10 side‐walls
};

// 7.4 Cylinder tops (at z = z₀)
Physical Surface("CylinderTop") = { 102, 103, 104, 105, 106, 107, 108, 109, 110, 111 };

// 7.5 Cylinder bottoms (at z = 0)
Physical Surface("CylinderBottom") = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

// 7.6 Channel top (ceiling) without holes
Physical Surface("ChannelTop") = { 101 };

// 7.7 Channel bottom (floor) without holes
Physical Surface("ChannelBottom") = { 1 };


// -----------------------------------------------------------------------------
// 8. Physical Points 
//-----------------------------------------------------------------------------

// 8.1 Inlet 
Physical Point("Inlet") = { 3, 4, 1003, 1004 };

// 8.2 Outlet 
Physical Point("Outlet") = {
  9, 10, 13, 14,         // bottom edge points
  1009, 1010, 1013, 1014 // corresponding top edge points
};

// 8.3 LateralWalls 
Physical Point("LateralWalls") = {
  1, 2, 5, 6, 7, 8, 11, 12, 15, 16,
  1001, 1002, 1005, 1006, 1007, 1008, 1011, 1012, 1015, 1016
};

// 8.4 Interface (all cylindrical side‐wall points, bottom and top circle‐arc endpoints)
Physical Point("Interf") = {
  // bottom arcs:
  101, 102, 103, 104, 106, 107, 108, 109,
  111, 112, 113, 114, 116, 117, 118, 119,
  121, 122, 123, 124, 126, 127, 128, 129,
  131, 132, 133, 134, 136, 137, 138, 139,
  141, 142, 143, 144, 146, 147, 148, 149,
  // top arcs:
  1101, 1102, 1103, 1104, 1106, 1107, 1108, 1109,
  1111, 1112, 1113, 1114, 1116, 1117, 1118, 1119,
  1121, 1122, 1123, 1124, 1126, 1127, 1128, 1129,
  1131, 1132, 1133, 1134, 1136, 1137, 1138, 1139,
  1141, 1142, 1143, 1144, 1146, 1147, 1148, 1149
};

// 8.5 CylinderTops (centers of each obstacle at z=z₀)
Physical Point("CylinderTop") = {
  1100, 1105, 1110, 1115, 1120,
  1125, 1130, 1135, 1140, 1145
};

// 8.6 CylinderBottoms (centers of each obstacle at z=0)
Physical Point("CylinderBottom") = {
  100, 105, 110, 115, 120, 125, 130, 135, 140, 145
};

// 8.7 ChannelTop (outer ceiling loop )
Physical Point("ChannelTop") = {
  1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008,
  1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016
};

// 8.8 ChannelBottom (outer floor loop )
Physical Point("ChannelBottom") = {
  1, 2, 3, 4, 5, 6, 7, 8,
  9, 10, 11, 12, 13, 14, 15, 16
};

// -----------------------------------------------------------------------------
// 9. Tag them as Physical Volumes
// -----------------------------------------------------------------------------

Physical Volume("PorousDomain") = { 5001, 5002, 5003, 5004, 5005, 5006, 5007, 5008, 5009, 5010 };
Physical Volume("FluidDomain")  = { 5100 };

// COMPLETED //
