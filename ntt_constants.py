# File generated with `python generate_constants.py`.
# Precomputations for NTT (reverse table of ψ and its inverse for different n.
ψ_rev = dict()
ψ = dict()
ψ_inv_rev = dict()
ψ_inv = dict()
n_inv = dict()
ψ_rev[12289] = [1, 1479, 4043, 7143, 5736, 4134, 1305, 722, 1646, 1212, 6429, 9094, 3504, 8747, 9744, 8668, 4591, 6561, 5023, 6461, 10938, 4978, 6512, 8961, 11340, 9664, 9650, 4821, 563, 9314, 2744, 3006, 1000, 4320, 12208, 3091, 9326, 4896, 2366, 9238, 11563, 7678, 1853, 140, 1635, 9521, 11112, 4255, 7203, 10963, 9088, 9275, 790, 955, 11119, 2319, 9542, 4846, 3135, 3712, 9995, 11227, 3553, 7484, 544, 5791, 11950, 2468, 11267, 9, 9447, 11809, 10616, 8011, 7300, 6958, 1381, 2525, 4177, 8705, 2837, 5374, 4354, 130, 2396, 4452, 3296, 8340, 12171, 9813, 2197, 5067, 11336, 3748, 5767, 827, 3284, 2881, 5092, 10200, 10276, 9000, 9048, 11560, 10593, 10861, 334, 2426, 4632, 5755, 11029, 4388, 10530, 3707, 3694, 7110, 11934, 3382, 2548, 8058, 4890, 6378, 9558, 3932, 5542, 12144, 3459, 3637, 1663, 1777, 1426, 7635, 2704, 5291, 7351, 8653, 9140, 160, 12286, 7852, 2166, 8374, 7370, 12176, 3364, 10600, 9018, 4057, 2174, 7917, 2847, 7875, 7094, 9509, 10805, 4895, 2305, 5042, 4053, 9644, 3985, 7384, 476, 3531, 420, 6730, 2178, 1544, 9273, 243, 9289, 11618, 3136, 5191, 8889, 9890, 9103, 6882, 10163, 1630, 11136, 2884, 8241, 10040, 3247, 9603, 2969, 3978, 6957, 3510, 9919, 9424, 7575, 8146, 1537, 12047, 8585, 2678, 5019, 545, 7404, 1017, 10657, 7205, 10849, 8526, 3066, 12262, 11244, 2859, 2481, 7277, 2912, 5698, 354, 7428, 390, 11516, 3778, 8456, 442, 2401, 5101, 11222, 4976, 10682, 875, 3780, 7278, 11287, 5088, 4284, 6022, 9302, 2437, 3646, 10102, 9723, 6039, 9867, 11854, 7952, 10911, 1912, 11796, 8193, 9908, 5444, 9041, 1207, 5277, 1168, 11885, 4645, 1065, 2143, 3957, 2839, 10162, 151, 11858, 1579, 2505, 5906, 52, 3174, 1323, 2766, 3336, 6055, 6415, 677, 3445, 7509, 4698, 5057, 12097, 10968, 10240, 4912, 5241, 9369, 3127, 4169, 3482, 787, 6821, 11279, 12231, 241, 11286, 3532, 11404, 6008, 10333, 7280, 2844, 3438, 8077, 975, 5681, 8812, 142, 1105, 4080, 421, 3602, 6221, 4624, 6212, 3263, 8689, 5886, 4782, 5594, 3029, 4213, 504, 605, 9987, 2033, 8291, 10367, 8410, 11316, 11035, 10930, 5435, 3710, 6196, 6950, 5446, 8301, 468, 11973, 11907, 6152, 4948, 11889, 10561, 6153, 6427, 3643, 5415, 56, 9090, 5206, 6760, 1702, 10302, 11635, 3565, 5315, 8214, 7373, 4324, 10120, 11767, 5079, 3262, 11011, 2344, 6715, 1973, 5925, 1018, 3514, 11248, 7500, 7822, 5537, 4749, 8500, 12142, 5456, 7840, 6844, 8429, 7753, 1050, 6118, 3818, 9606, 1190, 5876, 2281, 2031, 5333, 8298, 8320, 12133, 2767, 453, 6381, 418, 3772, 5429, 4774, 1293, 7552, 2361, 1843, 9259, 4115, 218, 2908, 8855, 8760, 2882, 10484, 1954, 2051, 2447, 6147, 576, 3963, 1858, 7535, 3315, 11863, 2925, 347, 3757, 1975, 10596, 3009, 174, 11566, 9551, 5868, 2655, 6554, 1512, 11939, 5383, 10474, 9087, 7796, 6920, 10232, 6374, 1483, 49, 11026, 1489, 2500, 10706, 5942, 1404, 11964, 11143, 948, 4049, 3728, 1159, 5990, 652, 5766, 6190, 11994, 4016, 4077, 2919, 3762, 6328, 7183, 10695, 1962, 7991, 8960, 12121, 9597, 7105, 1200, 6122, 9734, 3956, 1360, 6119, 5297, 3054, 6803, 9166, 1747, 5919, 4433, 3834, 5257, 683, 2459, 8633, 12225, 9786, 9341, 6507, 1566, 11454, 6224, 3570, 8049, 3150, 1319, 4046, 11580, 1958, 7967, 2078, 1112, 11231, 8210, 11367, 441, 1826,
                9363, 9118, 4489, 3708, 3238, 11153, 3449, 7080, 1092, 3359, 3205, 8024, 8611, 10361, 11825, 2068, 10900, 4404, 346, 3163, 8257, 7449, 6127, 12164, 11749, 10763, 4222, 8051, 11677, 8921, 8062, 7228, 11071, 11851, 3515, 9011, 5993, 6877, 8080, 1536, 10568, 4103, 9860, 11572, 8700, 1373, 2982, 3448, 11946, 4538, 1908, 4727, 11081, 1866, 7078, 10179, 716, 10125, 6873, 1705, 2450, 11475, 416, 10224, 5826, 7725, 8794, 1756, 4145, 8755, 8328, 5063, 4176, 8524, 10771, 2461, 2275, 8022, 5653, 6693, 6302, 11710, 3889, 212, 6323, 9175, 2769, 5734, 1176, 5508, 11014, 4860, 11164, 11158, 10844, 11841, 1014, 7508, 7365, 10962, 3607, 5232, 8347, 12221, 10029, 7723, 5836, 3200, 1535, 9572, 60, 7784, 10032, 10872, 5676, 3087, 6454, 7406, 3975, 7326, 8545, 2528, 3056, 5845, 5588, 11877, 5102, 1255, 506, 10897, 5784, 9615, 2212, 3338, 9013, 1178, 9513, 6811, 8778, 10347, 3408, 1165, 2575, 10453, 425, 11897, 10104, 377, 4578, 375, 1620, 1038, 11366, 6085, 4167, 6092, 2231, 2800, 12096, 1522, 2151, 8946, 8170, 5002, 12269, 7681, 5163, 10545, 1314, 2894, 3654, 11951, 3947, 9834, 6599, 7350, 7174, 1248, 2442, 8330, 6492, 6330, 10141, 5724, 10964, 1945, 1029, 8945, 6691, 10397, 3624, 6825, 4906, 4670, 512, 7735, 11295, 9389, 12050, 1804, 1403, 6195, 7100, 406, 10602, 7021, 12143, 8914, 9998, 7954, 3393, 8464, 8054, 7376, 8761, 11667, 1737, 4499, 5672, 8307, 9342, 11653, 5609, 4605, 2689, 180, 8151, 5219, 1409, 204, 6780, 9806, 2054, 1344, 9247, 463, 8882, 3981, 1468, 4475, 7043, 3017, 1236, 9168, 4705, 2600, 11232, 4739, 4251, 1226, 6771, 11925, 2360, 3028, 5216, 11839, 10345, 11711, 5368, 11779, 7628, 2622, 6903, 8929, 7605, 7154, 12226, 8481, 8619, 2373, 7302, 10891, 9199, 826, 5043, 5789, 8787, 6671, 10631, 9224, 1506, 7806, 5703, 4719, 11538, 6389, 11379, 4693, 9951, 11872, 9996, 6138, 8820, 4443, 8871, 7186, 10398, 1802, 10734, 1590, 4411, 1223, 2334, 2946, 6828, 2637, 4510, 881, 365, 10362, 1015, 7250, 6742, 2485, 904, 24, 10918, 11009, 11675, 980, 11607, 5082, 7699, 5207, 8239, 844, 7087, 3221, 8016, 8452, 2595, 5289, 6627, 567, 2941, 1406, 2633, 6940, 2945, 3232, 11996, 3769, 7434, 3944, 8190, 6759, 5604, 11024, 9282, 10118, 8809, 9169, 6184, 6643, 6086, 8753, 5370, 8348, 8536, 1282, 3572, 9457, 2021, 4730, 3229, 1706, 3929, 5054, 3154, 9004, 7929, 12282, 1936, 8566, 11444, 11520, 5526, 50, 216, 767, 3805, 4153, 10076, 1279, 11424, 9617, 5170, 12100, 3116, 10080, 1763, 3815, 1734, 1350, 5832, 8420, 4423, 1530, 1694, 10036, 10421, 9559, 5411, 4820, 1160, 9195, 7771, 2840, 9811, 4194, 9270, 7315, 4565, 7211, 10506, 944, 7519, 7002, 8620, 7624, 6883, 3020, 5673, 5410, 1251, 10499, 7014, 2035, 11249, 6164, 10407, 8176, 12217, 10447, 3840, 2712, 4834, 2828, 4352, 1241, 4378, 3451, 4094, 3045, 5781, 9646, 11194, 7592, 8711, 8823, 10588, 7785, 11511, 2626, 530, 10808, 9332, 9349, 2046, 8972, 9757, 8957, 12150, 3268, 3795, 1849, 6513, 4523, 4301, 457, 8, 8835, 3758, 8071, 4390, 10013, 982, 2593, 879, 9687, 10388, 11787, 7171, 6063, 8496, 8443, 1573, 5969, 4649, 9360, 6026, 1030, 11823, 10608, 8468, 11415, 9988, 5650, 12119, 648, 12139, 2307, 8000, 11498, 9855, 9416, 2827, 9754, 11169, 21, 6481]
ψ[12289] = [1, 1826, 3957, 11839, 1663, 1255, 5876, 1279, 544, 10224, 2033, 980, 7575, 6825, 1404, 7592, 1000, 7228, 12231, 4693, 3985, 1522, 1858, 944, 3284, 11841, 5315, 9169, 4976, 4605, 3054, 9687, 4591, 2068, 3445, 10891, 3364, 10453, 2361, 10036, 2837, 6693, 6152, 1406, 11244, 8914, 6328, 3268, 7203, 3448, 4080, 2946, 9103, 7350, 1512, 8176, 10530, 7784, 7500, 5054, 11854, 4475, 11454, 11415, 1646, 7080, 52, 8929, 9140, 1178, 453, 3815, 10616, 5063, 3710, 3221, 7404, 1804, 652, 10808, 11563, 1536, 2844, 7186, 9273, 10545, 10596, 5410, 10593, 12221, 11011, 1282, 6022, 9806, 683, 5969, 11340, 12164, 5241, 9224, 7094, 1038, 2882, 2840, 12171, 5734, 56, 3944, 390, 11667, 7105, 8835, 9542, 10179, 5886, 7250, 3247, 5724, 6374, 1241, 4890, 7326, 6844, 11520, 9041, 4739, 1958, 11498, 5736, 3708, 11858, 11779, 2704, 9615, 8298, 12100, 11267, 1756, 11316, 5207, 8585, 7735, 4049, 7785, 9326, 9011, 11404, 6138, 420, 5002, 2925, 7624, 10276, 10962, 10120, 8753, 7278, 5219, 5919, 6063, 10938, 3163, 12097, 5789, 2174, 377, 218, 4820, 2396, 212, 6153, 3232, 2912, 8464, 7991, 4523, 790, 4727, 4624, 881, 11136, 8330, 9087, 2712, 11934, 3087, 8500, 12282, 11796, 9168, 3150, 648, 3504, 8024, 3336, 8481, 2166, 10347, 5429, 8420, 1381, 2461, 8301, 5289, 10849, 406, 4016, 8972, 1635, 11572, 5681, 1590, 3136, 11951, 9551, 2035, 4632, 3200, 5925, 4730, 10102, 463, 9786, 1030, 563, 8051, 3482, 4719, 2305, 6092, 2447, 7315, 11336, 4860, 1702, 11024, 442, 8307, 3956, 10013, 9995, 1705, 4213, 24, 6957, 8945, 1489, 3045, 5542, 5845, 6118, 767, 11885, 11925, 11231, 9754, 4043, 9118, 10162, 11711, 1426, 10897, 2031, 9617, 11950, 7725, 10367, 5082, 1537, 4670, 11143, 8823, 12208, 11851, 11286, 11872, 476, 8946, 3315, 7002, 5092, 7508, 7373, 6643, 875, 180, 9166, 11787, 5023, 4404, 4698, 826, 9018, 11897, 9259, 9559, 4354, 11710, 11889, 6940, 2481, 7954, 10695, 1849, 9088, 4538, 3602, 2637, 10163, 1248, 5383, 10447, 3694, 10872, 5537, 9004, 10911, 3017, 3570, 5650, 6429, 3359, 1323, 7154, 12286, 6811, 418, 1350, 7300, 8524, 6950, 8452, 10657, 6195, 6190, 9349, 1853, 4103, 8077, 1802, 9289, 2894, 174, 10499, 334, 7723, 6715, 9457, 2437, 1344, 8633, 9360, 9650, 10763, 3127, 7806, 10805, 6085, 1954, 4194, 2197, 5508, 5206, 6759, 3778, 4499, 6122, 8071, 3135, 10125, 5594, 2485, 2969, 1945, 49, 3451, 9558, 2528, 7753, 50, 5277, 1226, 2078, 9416, 1305, 11153, 2505, 2622, 7351, 3338, 12133, 10080, 9447, 8755, 10930, 844, 5019, 9389, 1159, 2626, 2366, 6877, 10333, 4443, 2178, 7681, 3757, 3020, 9048, 5232, 5079, 8348, 5088, 204, 3834, 8443, 6512, 7449, 10240, 6671, 2847, 375, 8855, 9195, 3296, 9175, 3643, 3769, 354, 7376, 12121, 457, 11119, 1866, 3263, 10362, 8241, 6330, 6920, 2828, 2548, 7406, 5456, 8566, 9908, 2600, 4046, 2307, 9744, 10361, 6415, 2373, 7370, 1165, 1293, 1530, 4177, 8022, 11973, 567, 3066, 7021, 2919, 8957, 11112, 1373, 142, 1223, 8889, 9834, 2655, 6164, 11029, 9572, 3514, 1706, 6039, 3981, 6507, 10608, 2744, 8921, 6821, 6389, 4053, 2800, 576, 7211, 5767, 11158, 11635, 10118, 5101, 11653, 6119, 2593, 3553, 11475, 605, 11009, 9919, 10397, 10706, 9646, 3459, 11877, 9606, 4153, 1065, 3028, 11367,
            21, 1479, 9363, 2839, 10345, 1777, 506, 2281, 11424, 5791, 5826, 8291, 11607, 8146, 4906, 11964, 8711, 4320, 11071, 241, 9951, 7384, 2151, 7535, 7519, 2881, 1014, 8214, 6184, 10682, 2689, 6803, 10388, 6561, 10900, 7509, 9199, 10600, 425, 1843, 10421, 5374, 6302, 4948, 2633, 2859, 9998, 7183, 3795, 10963, 11946, 421, 6828, 6882, 7174, 11939, 12217, 3707, 10032, 7822, 3154, 7952, 7043, 6224, 9988, 1212, 1092, 3174, 7605, 160, 9513, 6381, 1734, 8011, 4176, 6196, 8016, 1017, 1403, 5766, 9332, 7678, 10568, 3438, 10398, 243, 1314, 3009, 1251, 10861, 10029, 2344, 3572, 9302, 2054, 2459, 4649, 9664, 11749, 9369, 1506, 9509, 11366, 10484, 9811, 9813, 1176, 9090, 8190, 11516, 1737, 1200, 3758, 4846, 716, 4782, 6742, 9603, 10964, 1483, 4378, 6378, 8545, 8429, 5526, 1207, 4251, 7967, 9855, 4134, 3238, 1579, 7628, 5291, 2212, 8320, 3116, 9, 4145, 11035, 8239, 2678, 11295, 3728, 11511, 4896, 5993, 6008, 8820, 6730, 12269, 347, 6883, 9000, 3607, 11767, 5370, 11287, 1409, 4433, 8496, 4978, 8257, 10968, 8787, 7917, 4578, 2908, 1160, 4452, 6323, 6427, 11996, 5698, 8054, 8960, 4301, 955, 11081, 6212, 365, 2884, 6492, 7796, 4834, 3382, 6454, 12142, 1936, 8193, 4705, 1319, 12139, 8747, 8611, 6055, 8619, 8374, 3408, 4774, 4423, 2525, 2275, 468, 6627, 8526, 10602, 4077, 9757, 9521, 8700, 8812, 4411, 5191, 3947, 5868, 11249, 5755, 1535, 1018, 3229, 9723, 8882, 9341, 11823, 9314, 11677, 787, 11538, 5042, 2231, 6147, 4565, 3748, 11164, 10302, 9282, 2401, 9342, 1360, 982, 11227, 2450, 504, 10918, 3510, 6691, 2500, 5781, 12144, 5588, 3818, 3805, 4645, 2360, 8210, 11169, 7143, 4489, 151, 5368, 7635, 5784, 5333, 5170, 2468, 8794, 8410, 7699, 12047, 512, 948, 10588, 3091, 3515, 3532, 9996, 3531, 8170, 11863, 8620, 10200, 7365, 4324, 6086, 3780, 8151, 1747, 7171, 6461, 346, 5057, 5043, 4057, 10104, 4115, 5411, 130, 3889, 10561, 2945, 7277, 3393, 1962, 6513, 9275, 1908, 6221, 4510, 1630, 2442, 10474, 3840, 7110, 5676, 4749, 7929, 1912, 1236, 8049, 12119, 9094, 3205, 2766, 12226, 7852, 8778, 3772, 5832, 6958, 10771, 5446, 2595, 7205, 7100, 11994, 2046, 140, 9860, 975, 10734, 11618, 3654, 11566, 7014, 2426, 5836, 1973, 2021, 3646, 9247, 12225, 6026, 4821, 4222, 4169, 5703, 4895, 4167, 2051, 9270, 5067, 11014, 6760, 5604, 8456, 5672, 9734, 4390, 3712, 6873, 3029, 904, 3978, 1029, 11026, 4094, 3932, 3056, 1050, 216, 1168, 6771, 1112, 2827, 722, 3449, 5906, 6903, 8653, 9013, 2767, 1763, 11809, 8328, 5435, 7087, 545, 12050, 5990, 530, 9238, 8080, 7280, 8871, 1544, 5163, 1975, 5673, 11560, 8347, 3262, 8536, 4284, 6780, 5257, 1573, 8961, 6127, 4912, 10631, 7875, 1620, 8760, 7771, 8340, 2769, 5415, 7434, 7428, 8761, 9597, 8, 2319, 7078, 8689, 1015, 10040, 10141, 10232, 4352, 8058, 3975, 7840, 11444, 5444, 11232, 11580, 8000, 8668, 11825, 677, 7302, 12176, 2575, 7552, 1694, 8705, 5653, 11907, 2941, 12262, 12143, 3762, 12150, 4255, 2982, 1105, 2334, 9890, 6599, 6554, 10407, 4388, 60, 11248, 3929, 9867, 1468, 1566, 8468, 3006, 8062, 11279, 11379, 9644, 12096, 3963, 10506, 827, 10844, 3565, 8809, 11222, 5609, 5297, 879, 7484, 416, 9987, 11675, 9424, 3624, 5942, 11194, 3637, 5102, 1190, 10076, 2143, 5216, 441, 6481]
ψ_inv_rev[12289] = [1, 10810, 5146, 8246, 11567, 10984, 8155, 6553, 3621, 2545, 3542, 8785, 3195, 5860, 11077, 10643, 9283, 9545, 2975, 11726, 7468, 2639, 2625, 949, 3328, 5777, 7311, 1351, 5828, 7266, 5728, 7698, 4805, 8736, 1062, 2294, 8577, 9154, 7443, 2747, 9970, 1170, 11334, 11499, 3014, 3201, 1326, 5086, 8034, 1177, 2768, 10654, 12149, 10436, 4611, 726, 3051, 9923, 7393, 2963, 9198, 81, 7969, 11289, 8652, 8830, 145, 6747, 8357, 2731, 5911, 7399, 4231, 9741, 8907, 355, 5179, 8595, 8582, 1759, 7901, 1260, 6534, 7657, 9863, 11955, 1428, 1696, 729, 3241, 3289, 2013, 2089, 7197, 9408, 9005, 11462, 6522, 8541, 953, 7222, 10092, 2476, 118, 3949, 8993, 7837, 9893, 12159, 7935, 6915, 9452, 3584, 8112, 9764, 10908, 5331, 4989, 4278, 1673, 480, 2842, 12280, 1022, 9821, 339, 6498, 11745, 10146, 11224, 7644, 404, 11121, 7012, 11082, 3248, 6845, 2381, 4096, 493, 10377, 1378, 4337, 435, 2422, 6250, 2566, 2187, 8643, 9852, 2987, 6267, 8005, 7201, 1002, 5011, 8509, 11414, 1607, 7313, 1067, 7188, 9888, 11847, 3833, 8511, 773, 11899, 4861, 11935, 6591, 9377, 5012, 9808, 9430, 1045, 27, 9223, 3763, 1440, 5084, 1632, 11272, 4885, 11744, 7270, 9611, 3704, 242, 10752, 4143, 4714, 2865, 2370, 8779, 5332, 8311, 9320, 2686, 9042, 2249, 4048, 9405, 1153, 10659, 2126, 5407, 3186, 2399, 3400, 7098, 9153, 671, 3000, 12046, 3016, 10745, 10111, 5559, 11869, 8758, 11813, 4905, 8304, 2645, 8236, 7247, 9984, 7394, 1484, 2780, 5195, 4414, 9442, 4372, 10115, 8232, 3271, 1689, 8925, 113, 4919, 3915, 10123, 4437, 3, 12129, 3149, 3636, 4938, 6998, 9585, 4654, 10863, 10512, 10626, 11848, 922, 4079, 1058, 11177, 10211, 4322, 10331, 709, 8243, 10970, 9139, 4240, 8719, 6065, 835, 10723, 5782, 2948, 2503, 64, 3656, 9830, 11606, 7032, 8455, 7856, 6370, 10542, 3123, 5486, 9235, 6992, 6170, 10929, 8333, 2555, 6167, 11089, 5184, 2692, 168, 3329, 4298, 10327, 1594, 5106, 5961, 8527, 9370, 8212, 8273, 295, 6099, 6523, 11637, 6299, 11130, 8561, 8240, 11341, 1146, 325, 10885, 6347, 1583, 9789, 10800, 1263, 12240, 10806, 5915, 2057, 5369, 4493, 3202, 1815, 6906, 350, 10777, 5735, 9634, 6421, 2738, 723, 12115, 9280, 1693, 10314, 8532, 11942, 9364, 426, 8974, 4754, 10431, 8326, 11713, 6142, 9842, 10238, 10335, 1805, 9407, 3529, 3434, 9381, 12071, 8174, 3030, 10446, 9928, 4737, 10996, 7515, 6860, 8517, 11871, 5908, 11836, 9522, 156, 3969, 3991, 6956, 10258, 10008, 6413, 11099, 2683, 8471, 6171, 11239, 4536, 3860, 5445, 4449, 6833, 147, 3789, 7540, 6752, 4467, 4789, 1041, 8775, 11271, 6364, 10316, 5574, 9945, 1278, 9027, 7210, 522, 2169, 7965, 4916, 4075, 6974, 8724, 654, 1987, 10587, 5529, 7083, 3199, 12233, 6874, 8646, 5862, 6136, 1728, 400, 7341, 6137, 382, 316, 11821, 3988, 6843, 5339, 6093, 8579, 6854, 1359, 1254, 973, 3879, 1922, 3998, 10256, 2302, 11684, 11785, 8076, 9260, 6695, 7507, 6403, 3600, 9026, 6077, 7665, 6068, 8687, 11868, 8209, 11184, 12147, 3477, 6608, 11314, 4212, 8851, 9445, 5009, 1956, 6281, 885, 8757, 1003, 12048, 58, 1010, 5468, 11502, 8807, 8120, 9162, 2920, 7048, 7377, 2049, 1321, 192, 7232, 7591, 4780, 8844, 11612, 5874, 6234, 8953, 9523, 10966, 9115, 12237, 6383, 9784, 10710, 431, 12138, 2127,
                    9450, 8332, 5808, 12268, 1120, 2535, 9462, 2873, 2434, 791, 4289, 9982, 150, 11641, 170, 6639, 2301, 874, 3821, 1681, 466, 11259, 6263, 2929, 7640, 6320, 10716, 3846, 3793, 6226, 5118, 502, 1901, 2602, 11410, 9696, 11307, 2276, 7899, 4218, 8531, 3454, 12281, 11832, 7988, 7766, 5776, 10440, 8494, 9021, 139, 3332, 2532, 3317, 10243, 2940, 2957, 1481, 11759, 9663, 778, 4504, 1701, 3466, 3578, 4697, 1095, 2643, 6508, 9244, 8195, 8838, 7911, 11048, 7937, 9461, 7455, 9577, 8449, 1842, 72, 4113, 1882, 6125, 1040, 10254, 5275, 1790, 11038, 6879, 6616, 9269, 5406, 4665, 3669, 5287, 4770, 11345, 1783, 5078, 7724, 4974, 3019, 8095, 2478, 9449, 4518, 3094, 11129, 7469, 6878, 2730, 1868, 2253, 10595, 10759, 7866, 3869, 6457, 10939, 10555, 8474, 10526, 2209, 9173, 189, 7119, 2672, 865, 11010, 2213, 8136, 8484, 11522, 12073, 12239, 6763, 769, 845, 3723, 10353, 7, 4360, 3285, 9135, 7235, 8360, 10583, 9060, 7559, 10268, 2832, 8717, 11007, 3753, 3941, 6919, 3536, 6203, 5646, 6105, 3120, 3480, 2171, 3007, 1265, 6685, 5530, 4099, 8345, 4855, 8520, 293, 9057, 9344, 5349, 9656, 10883, 9348, 11722, 5662, 7000, 9694, 3837, 4273, 9068, 5202, 11445, 4050, 7082, 4590, 7207, 682, 11309, 614, 1280, 1371, 12265, 11385, 9804, 5547, 5039, 11274, 1927, 11924, 11408, 7779, 9652, 5461, 9343, 9955, 11066, 7878, 10699, 1555, 10487, 1891, 5103, 3418, 7846, 3469, 6151, 2293, 417, 2338, 7596, 910, 5900, 751, 7570, 6586, 4483, 10783, 3065, 1658, 5618, 3502, 6500, 7246, 11463, 3090, 1398, 4987, 9916, 3670, 3808, 63, 5135, 4684, 3360, 5386, 9667, 4661, 510, 6921, 578, 1944, 450, 7073, 9261, 9929, 364, 5518, 11063, 8038, 7550, 1057, 9689, 7584, 3121, 11053, 9272, 5246, 7814, 10821, 8308, 3407, 11826, 3042, 10945, 10235, 2483, 5509, 12085, 10880, 7070, 4138, 12109, 9600, 7684, 6680, 636, 2947, 3982, 6617, 7790, 10552, 622, 3528, 4913, 4235, 3825, 8896, 4335, 2291, 3375, 146, 5268, 1687, 11883, 5189, 6094, 10886, 10485, 239, 2900, 994, 4554, 11777, 7619, 7383, 5464, 8665, 1892, 5598, 3344, 11260, 10344, 1325, 6565, 2148, 5959, 5797, 3959, 9847, 11041, 5115, 4939, 5690, 2455, 8342, 338, 8635, 9395, 10975, 1744, 7126, 4608, 20, 7287, 4119, 3343, 10138, 10767, 193, 9489, 10058, 6197, 8122, 6204, 923, 11251, 10669, 11914, 7711, 11912, 2185, 392, 11864, 1836, 9714, 11124, 8881, 1942, 3511, 5478, 2776, 11111, 3276, 8951, 10077, 2674, 6505, 1392, 11783, 11034, 7187, 412, 6701, 6444, 9233, 9761, 3744, 4963, 8314, 4883, 5835, 9202, 6613, 1417, 2257, 4505, 12229, 2717, 10754, 9089, 6453, 4566, 2260, 68, 3942, 7057, 8682, 1327, 4924, 4781, 11275, 448, 1445, 1131, 1125, 7429, 1275, 6781, 11113, 6555, 9520, 3114, 5966, 12077, 8400, 579, 5987, 5596, 6636, 4267, 10014, 9828, 1518, 3765, 8113, 7226, 3961, 3534, 8144, 10533, 3495, 4564, 6463, 2065, 11873, 814, 9839, 10584, 5416, 2164, 11573, 2110, 5211, 10423, 1208, 7562, 10381, 7751, 343, 8841, 9307, 10916, 3589, 717, 2429, 8186, 1721, 10753, 4209, 5412, 6296, 3278, 8774, 438, 1218, 5061, 4227, 3368, 612, 4238, 8067, 1526, 540, 125, 6162, 4840, 4032, 9126, 11943, 7885, 1389, 10221, 464, 1928, 3678, 4265, 9084, 8930, 11197, 5209, 8840, 1136, 9051, 8581, 7800, 3171, 2926, 10463]
ψ_inv[12289] = [1, 5808, 11848, 7073, 10146, 2213, 11099, 7187, 8652, 1095, 6347, 8665, 2865, 614, 2302, 11873, 4805, 11410, 6992, 6680, 1067, 3480, 8724, 1445, 11462, 1783, 8326, 193, 2645, 910, 1010, 4227, 9283, 3821, 10723, 10821, 2422, 8360, 1041, 12229, 7901, 1882, 5735, 5690, 2399, 9955, 11184, 9307, 8034, 139, 8527, 146, 27, 9348, 382, 6636, 3584, 10595, 4737, 9714, 113, 4987, 11612, 464, 3621, 4289, 709, 1057, 6845, 845, 4449, 8314, 4231, 7937, 2057, 2148, 2249, 11274, 3600, 5211, 9970, 12281, 2692, 3528, 4861, 4855, 6874, 9520, 3949, 4518, 3529, 10669, 4414, 1658, 7377, 6162, 3328, 10716, 7032, 5509, 8005, 3753, 9027, 3942, 729, 6616, 10314, 7126, 10745, 3418, 5009, 4209, 3051, 11759, 6299, 239, 11744, 5202, 6854, 3961, 480, 10526, 9522, 3276, 3636, 5386, 6383, 8840, 11567, 9462, 11177, 5518, 11121, 12073, 11239, 9233, 8357, 8195, 1263, 11260, 8311, 11385, 9260, 5416, 8577, 7899, 2555, 6617, 3833, 6685, 5529, 1275, 7222, 3019, 10238, 8122, 7394, 6586, 8120, 8067, 7468, 6263, 64, 3042, 8643, 10268, 10316, 6453, 9863, 5275, 723, 8635, 671, 1555, 11314, 2429, 12149, 10243, 295, 5189, 5084, 9694, 6843, 1518, 5331, 6457, 8517, 3511, 4437, 63, 9523, 9084, 3195, 170, 4240, 11053, 10377, 4360, 7540, 6613, 5179, 8449, 1815, 9847, 10659, 7779, 6068, 10381, 3014, 5776, 10327, 8896, 5012, 9344, 1728, 8400, 12159, 6878, 8174, 2185, 8232, 7246, 7232, 11943, 5828, 5118, 10542, 4138, 8509, 6203, 7965, 4924, 2089, 3669, 426, 4119, 8758, 2293, 8757, 8774, 9198, 1701, 11341, 11777, 242, 4590, 3879, 3495, 9821, 7119, 6956, 6505, 4654, 6921, 12138, 7800, 5146, 1120, 4079, 9929, 7644, 8484, 8471, 6701, 145, 6508, 9789, 5598, 8779, 1371, 11785, 9839, 1062, 11307, 10929, 2947, 9888, 3007, 1987, 1125, 8541, 7724, 6142, 10058, 7247, 751, 11502, 612, 2975, 466, 2948, 3407, 2566, 9060, 11271, 10754, 6534, 1040, 6421, 8342, 7098, 7878, 3477, 3589, 2768, 2532, 8212, 1687, 3763, 5662, 11821, 10014, 9764, 7866, 7515, 8881, 3915, 3670, 6234, 3678, 3542, 150, 10970, 7584, 4096, 10353, 147, 5835, 8907, 7455, 4493, 5797, 9405, 11924, 6077, 1208, 11334, 7988, 3329, 4235, 6591, 293, 5862, 5966, 7837, 11129, 9381, 7711, 4372, 3502, 1321, 4032, 7311, 3793, 7856, 10880, 1002, 6919, 522, 8682, 3289, 5406, 11942, 20, 5559, 3469, 6281, 6296, 7393, 778, 8561, 994, 9611, 4050, 1254, 8144, 12280, 9173, 3969, 10077, 6998, 4661, 10710, 9051, 8155, 2434, 4322, 8038, 11082, 6763, 3860, 3744, 5911, 7911, 10806, 1325, 2686, 5547, 7507, 11573, 7443, 8531, 11089, 10552, 773, 4099, 3199, 11113, 2476, 2478, 1805, 923, 2780, 10783, 2920, 540, 2625, 7640, 9830, 10235, 2987, 8717, 9945, 2260, 1428, 11038, 9280, 10975, 12046, 1891, 8851, 1721, 4611, 2957, 6523, 10886, 11272, 4273, 6093, 8113, 4278, 10555, 5908, 2776, 12129, 4684, 9115, 11197, 11077, 2301, 6065, 5246, 4337, 9135, 4467, 2257, 8582, 72, 350, 5115, 5407, 5461, 11868, 343, 1326, 8494, 5106, 2291, 9430, 9656, 7341, 5987, 6915, 1868, 10446, 11864, 1689, 3090, 4780, 1389, 5728, 1901, 5486, 9600, 1607, 6105, 4075, 11275, 9408, 4770, 4754, 10138, 4905, 2338, 12048, 1218, 7969, 3578, 325, 7383, 4143, 682, 3998, 6463, 6498, 865, 10008, 11783, 10512, 1944, 9450,
                2926, 10810, 12268, 922, 9261, 11224, 8136, 2683, 412, 8830, 2643, 1583, 1892, 2370, 1280, 11684, 814, 8736, 9696, 6170, 636, 7188, 2171, 654, 1131, 6522, 5078, 11713, 9489, 8236, 5900, 5468, 3368, 9545, 1681, 5782, 8308, 6250, 10583, 8775, 2717, 1260, 6125, 9634, 2455, 3400, 11066, 12147, 10916, 1177, 3332, 9370, 5268, 9223, 11722, 316, 4267, 8112, 10759, 10996, 11124, 4919, 9916, 5874, 1928, 2545, 9982, 8243, 9689, 2381, 3723, 6833, 4883, 9741, 9461, 5369, 5959, 4048, 1927, 9026, 10423, 1170, 11832, 168, 4913, 11935, 8520, 8646, 3114, 8993, 3094, 3434, 11914, 9442, 5618, 2049, 4840, 5777, 3846, 8455, 12085, 7201, 3941, 7210, 7057, 3241, 9269, 8532, 4608, 10111, 7846, 1956, 5412, 9923, 9663, 11130, 2900, 7270, 11445, 1359, 3534, 2842, 2209, 156, 8951, 4938, 9667, 9784, 1136, 10984, 2873, 10211, 11063, 7012, 12239, 4536, 9761, 2731, 8838, 12240, 10344, 9320, 9804, 6695, 2164, 9154, 4218, 6167, 7790, 8511, 5530, 7083, 6781, 10092, 8095, 10335, 6204, 1484, 4483, 9162, 1526, 2639, 2929, 3656, 10945, 9852, 2832, 5574, 4566, 11955, 1790, 12115, 9395, 3000, 10487, 4212, 8186, 10436, 2940, 6099, 6094, 1632, 3837, 5339, 3765, 4989, 10939, 11871, 5478, 3, 5135, 10966, 8930, 5860, 6639, 8719, 9272, 1378, 3285, 6752, 1417, 8595, 1842, 6906, 11041, 2126, 9652, 8687, 7751, 3201, 10440, 1594, 4335, 9808, 5349, 400, 579, 7935, 2730, 3030, 392, 3271, 11463, 7591, 7885, 7266, 502, 3123, 12109, 11414, 5646, 4916, 4781, 7197, 5287, 8974, 3343, 11813, 417, 1003, 438, 81, 3466, 1146, 7619, 10752, 7207, 1922, 4564, 339, 2672, 10258, 1392, 10863, 578, 2127, 3171, 8246, 2535, 1058, 364, 404, 11522, 6171, 6444, 6747, 9244, 10800, 3344, 5332, 12265, 8076, 10584, 2294, 2276, 8333, 3982, 11847, 1265, 10587, 7429, 953, 4974, 9842, 6197, 9984, 7570, 8807, 4238, 11726, 11259, 2503, 11826, 2187, 7559, 6364, 9089, 7657, 10254, 2738, 338, 9153, 10699, 6608, 717, 10654, 3317, 8273, 11883, 1440, 7000, 3988, 9828, 10908, 3869, 6860, 1942, 10123, 3808, 8953, 4265, 8785, 11641, 9139, 3121, 493, 7, 3789, 9202, 355, 9577, 3202, 3959, 1153, 11408, 7665, 7562, 11499, 7766, 4298, 3825, 9377, 9057, 6136, 12077, 9893, 7469, 12071, 11912, 10115, 6500, 192, 9126, 1351, 6226, 6370, 7070, 5011, 3536, 2169, 1327, 2013, 4665, 9364, 7287, 11869, 6151, 885, 3278, 2963, 4504, 8240, 4554, 3704, 7082, 973, 10533, 1022, 189, 3991, 2674, 9585, 510, 431, 8581, 6553, 791, 10331, 7550, 3248, 769, 5445, 4963, 7399, 11048, 5915, 6565, 9042, 5039, 6403, 2110, 2747, 3454, 5184, 622, 11899, 8345, 12233, 6555, 118, 9449, 9407, 11251, 5195, 3065, 7048, 125, 949, 6320, 11606, 2483, 6267, 11007, 1278, 68, 1696, 6879, 1693, 1744, 3016, 5103, 9445, 10753, 726, 1481, 11637, 10485, 4885, 9068, 8579, 7226, 1673, 8474, 11836, 11111, 3149, 3360, 12237, 5209, 10643, 874, 835, 7814, 435, 7235, 4789, 4505, 1759, 4113, 10777, 4939, 3186, 9343, 8209, 8841, 5086, 9021, 5961, 3375, 1045, 10883, 6137, 5596, 9452, 2253, 9928, 1836, 8925, 1398, 8844, 10221, 7698, 2602, 9235, 7684, 7313, 3120, 6974, 448, 9005, 11345, 10431, 10767, 8304, 7596, 58, 5061, 11289, 4697, 10885, 5464, 4714, 11309, 10256, 2065, 11745, 11010, 6413, 11034, 10626, 450, 8332, 10463]
# inverse of powers of 2 mod q
n_inv[12289] = {
    2: 6145,
    4: 9217,
    8: 10753,
    16: 11521,
    32: 11905,
    64: 12097,
    128: 12193,
    256: 12241,
    512: 12265,
    1024: 12277,
}
ψ_rev[3329] = [1, 1729, 2580, 3289, 2642, 630, 1897, 848, 1062, 1919, 193, 797, 2786, 3260, 569, 1746, 296, 2447, 1339, 1476, 3046, 56, 2240, 1333, 1426, 2094, 535, 2882, 2393, 2879, 1974, 821, 289, 331, 3253, 1756, 1197, 2304, 2277, 2055, 650, 1977, 2513, 632, 2865, 33, 1320, 1915, 2319, 1435, 807, 452, 1438, 2868, 1534, 2402, 2647, 2617, 1481, 648, 2474, 3110, 1227,
               910, 17, 2761, 583, 2649, 1637, 723, 2288, 1100, 1409, 2662, 3281, 233, 756, 2156, 3015, 3050, 1703, 1651, 2789, 1789, 1847, 952, 1461, 2687, 939, 2308, 2437, 2388, 733, 2337, 268, 641, 1584, 2298, 2037, 3220, 375, 2549, 2090, 1645, 1063, 319, 2773, 757, 2099, 561, 2466, 2594, 2804, 1092, 403, 1026, 1143, 2150, 2775, 886, 1722, 1212, 1874, 1029, 2110, 2935, 885, 2154]
ψ[3329] = [1, 17, 289, 1584, 296, 1703, 2319, 2804, 1062, 1409, 650, 1063, 1426, 939, 2647, 1722, 2642, 1637, 1197, 375, 3046, 1847, 1438, 1143, 2786, 756, 2865, 2099, 2393, 733, 2474, 2110, 2580, 583, 3253, 2037, 1339, 2789, 807, 403, 193, 3281, 2513, 2773, 535, 2437, 1481, 1874, 1897, 2288, 2277, 2090, 2240, 1461, 1534, 2775, 569, 3015, 1320, 2466, 1974, 268, 1227,
           885, 1729, 2761, 331, 2298, 2447, 1651, 1435, 1092, 1919, 2662, 1977, 319, 2094, 2308, 2617, 1212, 630, 723, 2304, 2549, 56, 952, 2868, 2150, 3260, 2156, 33, 561, 2879, 2337, 3110, 2935, 3289, 2649, 1756, 3220, 1476, 1789, 452, 1026, 797, 233, 632, 757, 2882, 2388, 648, 1029, 848, 1100, 2055, 1645, 1333, 2687, 2402, 886, 1746, 3050, 1915, 2594, 821, 641, 910, 2154]
ψ_inv_rev[3329] = [1, 1600, 40, 749, 2481, 1432, 2699, 687, 1583, 2760, 69, 543, 2532, 3136, 1410, 2267, 2508, 1355, 450, 936, 447, 2794, 1235, 1903, 1996, 1089, 3273, 283, 1853, 1990, 882, 3033, 2419, 2102, 219, 855, 2681, 1848, 712, 682, 927, 1795, 461, 1891, 2877, 2522, 1894, 1010, 1414, 2009, 3296, 464, 2697, 816, 1352, 2679, 1274, 1052, 1025, 2132, 1573, 76, 2998,
                   3040, 1175, 2444, 394, 1219, 2300, 1455, 2117, 1607, 2443, 554, 1179, 2186, 2303, 2926, 2237, 525, 735, 863, 2768, 1230, 2572, 556, 3010, 2266, 1684, 1239, 780, 2954, 109, 1292, 1031, 1745, 2688, 3061, 992, 2596, 941, 892, 1021, 2390, 642, 1868, 2377, 1482, 1540, 540, 1678, 1626, 279, 314, 1173, 2573, 3096, 48, 667, 1920, 2229, 1041, 2606, 1692, 680, 2746, 568, 3312]
ψ_inv[3329] = [1, 1175, 2419, 2688, 2508, 735, 1414, 279, 1583, 2443, 927, 642, 1996, 1684, 1274, 2229, 2481, 2300, 2681, 941, 447, 2572, 2697, 3096, 2532, 2303, 2877, 1540, 1853, 109, 1573, 680, 40, 394, 219, 992, 450, 2768, 3296, 1173, 69, 1179, 461, 2377, 3273, 780, 1025, 2606, 2699, 2117, 712, 1021, 1235, 3010, 1352, 667, 1410, 2237, 1894, 1678, 882, 1031, 2998,
               568, 1600, 2444, 2102, 3061, 1355, 863, 2009, 314, 2760, 554, 1795, 1868, 1089, 1239, 1052, 1041, 1432, 1455, 1848, 892, 2794, 556, 816, 48, 3136, 2926, 2522, 540, 1990, 1292, 76, 2746, 749, 1219, 855, 2596, 936, 1230, 464, 2573, 543, 2186, 1891, 1482, 283, 2954, 2132, 1692, 687, 1607, 682, 2390, 1903, 2266, 2679, 1920, 2267, 525, 1010, 1626, 3033, 1745, 3040, 3312]
# inverse of powers of 2 mod q
n_inv[3329] = {
    2: 1665,
    4: 2497,
    8: 2913,
    16: 3121,
    32: 3225,
    64: 3277,
    128: 3303,
    256: 3316,
    512: 1658,
    1024: 829,
}
