\ENCODING=ISO-8859-1
\Problem name: MCF_Model

Minimize
 obj: 30 F_01 + 34 F_02 + 68 F_03 + 69 F_04 + 25 F_05 + 36 F_06 + 27 F_07
      + 86 F_08 + 69 F_09 + 64 F_10 + 68 F_11 + 25 F_12 + 44 F_13 + 21 F_14
      + 34 F_15 + 13 F_16 + 42 F_17 + 19 F_18 + 39 F_19 + 76 F_20 + 24 F_21
      + 52 F_22 + 34 F_23 + 30 F_24 + 50 F_25 + 168 F_26 + 96 F_27 + 104 F_28
      + 110 F_29 + 30 F_30 + 93 F_31 + 30 F_32 + 152 F_33 + 69 F_34 + 100 F_35
      + 35 F_36 + 63 F_37 + 88 F_38 + 195 F_39 + 84 F_40 + 10000 S_01
      + 10000 S_02 + 10000 S_03 + 10000 S_04 + 10000 S_05 + 10000 S_06
      + 10000 S_07 + 10000 S_08 + 10000 S_09 + 10000 S_10 + 10000 S_11
      + 10000 S_12 + 10000 S_13 + 10000 S_14 + 10000 S_15 + 10000 S_16
      + 10000 S_17 + 10000 S_18 + 10000 S_19 + 10000 S_20 + 10000 S_21
      + 10000 S_22 + 10000 S_23 + 10000 S_24 + 10000 S_25 + 10000 S_26
      + 10000 S_27 + 10000 S_28 + 10000 S_29 + 10000 S_30 + 10000 S_31
      + 10000 S_32 + 10000 S_33 + 10000 S_34 + 10000 S_35 + 10000 S_36
      + 10000 S_37 + 10000 S_38 + 10000 S_39 + 10000 S_40 + 10000 S_41
      + 10000 S_42 + 10000 S_43 + 10000 S_44 + 10000 S_45 + 10000 S_46
      + 10000 S_47 + 10000 S_48 + 10000 S_49 + 10000 S_50 + 10000 S_51
      + 10000 S_52 + 10000 S_53 + 10000 S_54 + 10000 S_55 + 10000 S_56
      + 10000 S_57 + 10000 S_58 + 10000 S_59 + 10000 S_60
Subject To
 Demand_Verification_1:  F_01  = 1
 Demand_Verification_2:  F_02  = 1
 Demand_Verification_3:  F_03  = 1
 Demand_Verification_4:  F_04  = 1
 Demand_Verification_5:  F_05  = 1
 Demand_Verification_6:  F_06  = 1
 Demand_Verification_7:  F_07  = 1
 Demand_Verification_8:  F_08  = 1
 Demand_Verification_9:  F_09  = 1
 Demand_Verification_10: F_10  = 1
 Demand_Verification_11: F_11  = 1
 Demand_Verification_12: F_12  = 1
 Demand_Verification_13: F_13  = 1
 Demand_Verification_14: F_14  = 1
 Demand_Verification_15: F_15  = 1
 Demand_Verification_16: F_16  = 1
 Demand_Verification_17: F_17  = 1
 Demand_Verification_18: F_18  = 1
 Demand_Verification_19: F_19  = 1
 Demand_Verification_20: F_20  = 1
 Demand_Verification_21: F_21  = 1
 Demand_Verification_22: F_22  = 1
 Demand_Verification_23: F_23  = 1
 Demand_Verification_24: F_24  = 1
 Demand_Verification_25: F_25  = 1
 Demand_Verification_26: F_26  = 1
 Demand_Verification_27: F_27  = 1
 Demand_Verification_28: F_28  = 1
 Demand_Verification_29: F_29  = 1
 Demand_Verification_30: F_30  = 1
 Demand_Verification_31: F_31  = 1
 Demand_Verification_32: F_32  = 1
 Demand_Verification_33: F_33  = 1
 Demand_Verification_34: F_34  = 1
 Demand_Verification_35: F_35  = 1
 Demand_Verification_36: F_36  = 1
 Demand_Verification_37: F_37  = 1
 Demand_Verification_38: F_38  = 1
 Demand_Verification_39: F_39  = 1
 Demand_Verification_40: F_40  = 1
 Capacity_Constraint_1:  F_05  = 0
 Capacity_Constraint_2:  2 F_02 - RgCapacity_Constraint_2  = 0
 Capacity_Constraint_3:  - RgCapacity_Constraint_3  = 0
 Capacity_Constraint_4:  3 F_04 - RgCapacity_Constraint_4  = 0
 Capacity_Constraint_5:  2 F_06 - RgCapacity_Constraint_5  = 0
 Capacity_Constraint_6:  3 F_01 + 2 F_02 + 3 F_04 + F_12
                         - RgCapacity_Constraint_6  = 0
 Capacity_Constraint_7:   = 0
 Capacity_Constraint_8:  4 F_03  = 0
 Capacity_Constraint_9:  F_05  = 0
 Capacity_Constraint_10:  = 0
 Capacity_Constraint_11: F_07 + 2 F_08  = 0
 Capacity_Constraint_12: 2 F_32  = 0
 Capacity_Constraint_13: 2 F_20 + 4 F_26 + 5 F_39  = 0
 Capacity_Constraint_14:  = 0
 Capacity_Constraint_15: F_07 + 2 F_30  = 0
 Capacity_Constraint_16: 2 F_08 + 5 F_29 + 3 F_31  = 0
 Capacity_Constraint_17: - RgCapacity_Constraint_17  = 0
 Capacity_Constraint_18: 4 F_26 + 5 F_39  = 0
 Capacity_Constraint_19: 4 F_33 + 3 F_34  = 0
 Capacity_Constraint_20:  = 0
 Capacity_Constraint_21: 5 F_35  = 0
 Capacity_Constraint_22: 7 F_36 - RgCapacity_Constraint_22  = 0
 Capacity_Constraint_23: - RgCapacity_Constraint_23  = 0
 Capacity_Constraint_24: - RgCapacity_Constraint_24  = 0
 Capacity_Constraint_25: F_07 + 2 F_30 + 3 F_37  = 0
 Capacity_Constraint_26: 2 F_15 + 2 F_22 + 2 F_25  = 0
 Capacity_Constraint_27: 2 F_20  = 0
 Capacity_Constraint_28:  = 0
 Capacity_Constraint_29: 3 F_37  = 0
 Capacity_Constraint_30: 3 F_27  = 0
 Capacity_Constraint_31: 2 F_17 + F_18  = 0
 Capacity_Constraint_32: 3 F_09  = 0
 Capacity_Constraint_33: 2 F_08 + 3 F_19 + 3 F_31 + 3 F_40
                         - RgCapacity_Constraint_33  = 0
 Capacity_Constraint_34: 2 F_10 + F_16  = 0
 Capacity_Constraint_35: 2 F_23  = 0
 Capacity_Constraint_36: F_12 - RgCapacity_Constraint_36  = 0
 Capacity_Constraint_37:  = 0
 Capacity_Constraint_38: 4 F_11 + 4 F_33  = 0
 Capacity_Constraint_39: 2 F_13 + F_14  = 0
 Capacity_Constraint_40: 3 F_09 + 2 F_10  = 0
 Capacity_Constraint_41:  = 0
 Capacity_Constraint_42: F_14 + 2 F_17 + 4 F_26 + 4 F_28 + 5 F_39  = 0
 Capacity_Constraint_43:  = 0
 Capacity_Constraint_44:  = 0
 Capacity_Constraint_45: 3 F_24 + 4 F_33 + 4 F_38  = 0
 Capacity_Constraint_46: 3 F_34  = 0
 Capacity_Constraint_47: 2 F_13 + F_21  = 0
 Capacity_Constraint_48: F_18  = 0
 Capacity_Constraint_49: 3 F_19 + 2 F_20 - RgCapacity_Constraint_49  = 0
 Capacity_Constraint_50: F_21  = 0
 Capacity_Constraint_51: 2 F_22 + 2 F_23 - RgCapacity_Constraint_51  = 0
 Capacity_Constraint_52: 3 F_27 + 3 F_40 - 3 S_52  = 0
 Capacity_Constraint_53: 4 F_26 + 4 F_28 - RgCapacity_Constraint_53  = 0
 Capacity_Constraint_54: - RgCapacity_Constraint_54  = 0
 Capacity_Constraint_55: 2 F_25 + 3 F_37 - RgCapacity_Constraint_55  = 0
 Capacity_Constraint_56: 4 F_38 + 5 F_39  = 0
 Capacity_Constraint_57:  = 0
 Capacity_Constraint_58: 4 F_26 + 4 F_28  = 0
 Capacity_Constraint_59: 3 F_40 - 3 S_59  = 0
 Capacity_Constraint_60:  = 0
Bounds
      F_01 >= 0
      F_02 >= 0
      F_03 >= 0
      F_04 >= 0
      F_05 >= 0
      F_06 >= 0
      F_07 >= 0
      F_08 >= 0
      F_09 >= 0
      F_10 >= 0
      F_11 >= 0
      F_12 >= 0
      F_13 >= 0
      F_14 >= 0
      F_15 >= 0
      F_16 >= 0
      F_17 >= 0
      F_18 >= 0
      F_19 >= 0
      F_20 >= 0
      F_21 >= 0
      F_22 >= 0
      F_23 >= 0
      F_24 >= 0
      F_25 >= 0
      F_26 >= 0
      F_27 >= 0
      F_28 >= 0
      F_29 >= 0
      F_30 >= 0
      F_31 >= 0
      F_32 >= 0
      F_33 >= 0
      F_34 >= 0
      F_35 >= 0
      F_36 >= 0
      F_37 >= 0
      F_38 >= 0
      F_39 >= 0
      F_40 >= 0
      S_01 >= 0
      S_02 >= 0
      S_03 >= 0
      S_04 >= 0
      S_05 >= 0
      S_06 >= 0
      S_07 >= 0
      S_08 >= 0
      S_09 >= 0
      S_10 >= 0
      S_11 >= 0
      S_12 >= 0
      S_13 >= 0
      S_14 >= 0
      S_15 >= 0
      S_16 >= 0
      S_17 >= 0
      S_18 >= 0
      S_19 >= 0
      S_20 >= 0
      S_21 >= 0
      S_22 >= 0
      S_23 >= 0
      S_24 >= 0
      S_25 >= 0
      S_26 >= 0
      S_27 >= 0
      S_28 >= 0
      S_29 >= 0
      S_30 >= 0
      S_31 >= 0
      S_32 >= 0
      S_33 >= 0
      S_34 >= 0
      S_35 >= 0
      S_36 >= 0
      S_37 >= 0
      S_38 >= 0
      S_39 >= 0
      S_40 >= 0
      S_41 >= 0
      S_42 >= 0
      S_43 >= 0
      S_44 >= 0
      S_45 >= 0
      S_46 >= 0
      S_47 >= 0
      S_48 >= 0
      S_49 >= 0
      S_50 >= 0
      S_51 >= 0
      S_52 >= 0
      S_53 >= 0
      S_54 >= 0
      S_55 >= 0
      S_56 >= 0
      S_57 >= 0
      S_58 >= 0
      S_59 >= 0
      S_60 >= 0
 0 <= RgCapacity_Constraint_2 <= 5
 0 <= RgCapacity_Constraint_3 <= 5
 0 <= RgCapacity_Constraint_4 <= 2
 0 <= RgCapacity_Constraint_5 <= 3
 0 <= RgCapacity_Constraint_6 <= 4
 0 <= RgCapacity_Constraint_17 <= 5
 0 <= RgCapacity_Constraint_22 <= 4
 0 <= RgCapacity_Constraint_23 <= 2
 0 <= RgCapacity_Constraint_24 <= 2
 0 <= RgCapacity_Constraint_33 <= 5
 0 <= RgCapacity_Constraint_36 <= 8
 0 <= RgCapacity_Constraint_49 <= 2
 0 <= RgCapacity_Constraint_51 <= 8
 0 <= RgCapacity_Constraint_53 <= 5
 0 <= RgCapacity_Constraint_54 <= 10
 0 <= RgCapacity_Constraint_55 <= 5
Generals
 F_01  F_02  F_03  F_04  F_05  F_06  F_07  F_08  F_09  F_10  F_11  F_12  F_13 
 F_14  F_15  F_16  F_17  F_18  F_19  F_20  F_21  F_22  F_23  F_24  F_25  F_26 
 F_27  F_28  F_29  F_30  F_31  F_32  F_33  F_34  F_35  F_36  F_37  F_38  F_39 
 F_40  S_01  S_02  S_03  S_04  S_05  S_06  S_07  S_08  S_09  S_10  S_11  S_12 
 S_13  S_14  S_15  S_16  S_17  S_18  S_19  S_20  S_21  S_22  S_23  S_24  S_25 
 S_26  S_27  S_28  S_29  S_30  S_31  S_32  S_33  S_34  S_35  S_36  S_37  S_38 
 S_39  S_40  S_41  S_42  S_43  S_44  S_45  S_46  S_47  S_48  S_49  S_50  S_51 
 S_52  S_53  S_54  S_55  S_56  S_57  S_58  S_59  S_60 
End
