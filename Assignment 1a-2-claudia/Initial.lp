\ENCODING=ISO-8859-1
\Problem name: Initial

Minimize
 obj: 75 x1 + 30 x2 + 70 x3 + 25 x4 + 1000 x5 + 1000 x6 + 1000 x7 + 1000 x8
      + 1000 x9 + 1000 x10 + 1000 x11
Subject To
 Commodity_1: x1  = 1
 c2:          x2  = 1
 c3:          x3  = 1
 c4:          x4  = 1
 c5:          15 x1 - x5 - Rgc5  = 0
 c6:          5 x2 - x6 - Rgc6  = 0
 c7:          10 x3 - x7 - Rgc7  = 0
 c8:          15 x1 - x8 - Rgc8  = 0
 c9:          - x9 - Rgc9  = 0
 c10:         5 x2 + 10 x3 + 5 x4 - x10 - Rgc10  = 0
 c11:         - x11 - Rgc11  = 0
Bounds
 0 <= Rgc5 <= 20
 0 <= Rgc6 <= 10
 0 <= Rgc7 <= 10
 0 <= Rgc8 <= 20
 0 <= Rgc9 <= 40
 0 <= Rgc10 <= 10
 0 <= Rgc11 <= 30
End
