\ENCODING=ISO-8859-1
\Problem name: Initial

Minimize
 obj: 75 x1 + 30 x2 + 40 x3 + 70 x4 + 70 x5 + 25 x6 + 55 x7 + 1000 x8 + 1000 x9
      + 1000 x10 + 1000 x11 + 1000 x12 + 1000 x13 + 1000 x14
Subject To
 Commodity_1: x1  = 1
 c2:          x2 + x3  = 1
 c3:          x4 + x5  = 1
 c4:          x6 + x7  = 1
 c5:          15 x1 + 5 x3 - x8 - Rgc5  = 0
 c6:          5 x2 - x9 - Rgc6  = 0
 c7:          10 x4 - x10 - Rgc7  = 0
 c8:          15 x1 + 5 x3 + 10 x5 - x11 - Rgc8  = 0
 c9:          5 x7 - x12 - Rgc9  = 0
 c10:         5 x2 + 10 x4 + 5 x6 - x13 - Rgc10  = 0
 c11:         5 x3 + 10 x5 + 5 x7 - x14 - Rgc11  = 0
Bounds
 0 <= Rgc5 <= 20
 0 <= Rgc6 <= 10
 0 <= Rgc7 <= 10
 0 <= Rgc8 <= 20
 0 <= Rgc9 <= 40
 0 <= Rgc10 <= 10
 0 <= Rgc11 <= 30
End
