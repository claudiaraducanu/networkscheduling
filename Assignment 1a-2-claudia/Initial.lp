\ENCODING=ISO-8859-1
\Problem name: Initial

Minimize
 obj: 75 F_01 + 30 F_02 + 70 F_03 + 25 F_04 + 1000 S_01 + 1000 S_02 + 1000 S_03
      + 1000 S_04 + 1000 S_05 + 1000 S_06 + 1000 S_07
Subject To
 Bundle_1:              F_01  = 1
 Bundle_2:              F_02  = 1
 Bundle_3:              F_03  = 1
 Bundle_4:              F_04  = 1
 Capacity_Constraint_1: 15 F_01 - S_01 - RgCapacity_Constraint_1  = 0
 Capacity_Constraint_2: 5 F_02 - S_02 - RgCapacity_Constraint_2  = 0
 Capacity_Constraint_3: 10 F_03 - S_03 - RgCapacity_Constraint_3  = 0
 Capacity_Constraint_4: 15 F_01 - S_04 - RgCapacity_Constraint_4  = 0
 Capacity_Constraint_5: - S_05 - RgCapacity_Constraint_5  = 0
 Capacity_Constraint_6: 5 F_02 + 10 F_03 + 5 F_04 - S_06
                        - RgCapacity_Constraint_6  = 0
 Capacity_Constraint_7: - S_07 - RgCapacity_Constraint_7  = 0
Bounds
      F_01 >= 0
      F_02 >= 0
      F_03 >= 0
      F_04 >= 0
      S_01 >= 0
      S_02 >= 0
      S_03 >= 0
      S_04 >= 0
      S_05 >= 0
      S_06 >= 0
      S_07 >= 0
 0 <= RgCapacity_Constraint_1 <= 20
 0 <= RgCapacity_Constraint_2 <= 10
 0 <= RgCapacity_Constraint_3 <= 10
 0 <= RgCapacity_Constraint_4 <= 20
 0 <= RgCapacity_Constraint_5 <= 40
 0 <= RgCapacity_Constraint_6 <= 10
 0 <= RgCapacity_Constraint_7 <= 30
Generals
 F_01  F_02  F_03  F_04  S_01  S_02  S_03  S_04  S_05  S_06  S_07 
End
