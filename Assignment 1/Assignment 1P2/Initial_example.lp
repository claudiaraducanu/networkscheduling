\ENCODING=ISO-8859-1
\Problem name: Initial_example

Minimize
 obj: 150 x1 + 130 x2 + 150 x3 + 130 x4 + 150 x5 + 150 x6 + 100 x7 + 100 x8
      + 110 x9
Subject To
 Capacity_1: x1 + x2 + x9 >= 22
 c2:         x3 + x4 >= -4
 c3:         x2 + x5 + x9 >= -4
 c4:         x4 + x6 >= 22
 c5:         x7 >= 12
 c6:         x8 - 0.2 x9 >= -8
End
