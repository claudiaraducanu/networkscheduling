\ENCODING=ISO-8859-1
\Problem name: Initial

Minimize
 obj: 75 x1 + 30 x2 + 40 x3 + 70 x4 + 70 x5 + 25 x6 + 55 x7 + 1000 x8
Subject To
 c1:  x1  = 1
 c2:  x2 + x3  = 1
 c3:  x4 + x5  = 1
 c4:  x6 + x7  = 1
 c5:  15 x1 + 5 x3 - x9  = 0
 c6:  5 x2 - x10  = 0
 c7:  10 x4 - x11  = 0
 c8:  15 x1 + 5 x3 + 10 x5 - x12  = 0
 c9:  5 x7 - x13  = 0
 c10: 5 x2 + 10 x4 + 5 x6 - x8 - x14  = 0
 c11: 5 x3 + 10 x5 + 5 x7 - x15  = 0
Bounds
 0 <= x9 <= 20
 0 <= x10 <= 10
 0 <= x11 <= 10
 0 <= x12 <= 20
 0 <= x13 <= 40
 0 <= x14 <= 10
 0 <= x15 <= 30
End
