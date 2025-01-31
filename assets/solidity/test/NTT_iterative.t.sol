// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.25;

import {Test, console} from "forge-std/Test.sol";
import  "../src/NTT_Iterative.sol";



contract NTTTest is Test {
   
    NTT ntt=new NTT();
     uint256[1024] psi_rev = [uint256(1), 1479, 4043, 7143, 5736, 4134, 1305, 722, 1646, 1212, 6429, 9094, 3504, 8747, 9744, 8668, 4591, 6561, 5023, 6461, 10938, 4978, 6512, 8961, 11340, 9664, 9650, 4821, 563, 9314, 2744, 3006, 1000, 4320, 12208, 3091, 9326, 4896, 2366, 9238, 11563, 7678, 1853, 140, 1635, 9521, 11112, 4255, 7203, 10963, 9088, 9275, 790, 955, 11119, 2319, 9542, 4846, 3135, 3712, 9995, 11227, 3553, 7484, 544, 5791, 11950, 2468, 11267, 9, 9447, 11809, 10616, 8011, 7300, 6958, 1381, 2525, 4177, 8705, 2837, 5374, 4354, 130, 2396, 4452, 3296, 8340, 12171, 9813, 2197, 5067, 11336, 3748, 5767, 827, 3284, 2881, 5092, 10200, 10276, 9000, 9048, 11560, 10593, 10861, 334, 2426, 4632, 5755, 11029, 4388, 10530, 3707, 3694, 7110, 11934, 3382, 2548, 8058, 4890, 6378, 9558, 3932, 5542, 12144, 3459, 3637, 1663, 1777, 1426, 7635, 2704, 5291, 7351, 8653, 9140, 160, 12286, 7852, 2166, 8374, 7370, 12176, 3364, 10600, 9018, 4057, 2174, 7917, 2847, 7875, 7094, 9509, 10805, 4895, 2305, 5042, 4053, 9644, 3985, 7384, 476, 3531, 420, 6730, 2178, 1544, 9273, 243, 9289, 11618, 3136, 5191, 8889, 9890, 9103, 6882, 10163, 1630, 11136, 2884, 8241, 10040, 3247, 9603, 2969, 3978, 6957, 3510, 9919, 9424, 7575, 8146, 1537, 12047, 8585, 2678, 5019, 545, 7404, 1017, 10657, 7205, 10849, 8526, 3066, 12262, 11244, 2859, 2481, 7277, 2912, 5698, 354, 7428, 390, 11516, 3778, 8456, 442, 2401, 5101, 11222, 4976, 10682, 875, 3780, 7278, 11287, 5088, 4284, 6022, 9302, 2437, 3646, 10102, 9723, 6039, 9867, 11854, 7952, 10911, 1912, 11796, 8193, 9908, 5444, 9041, 1207, 5277, 1168, 11885, 4645, 1065, 2143, 3957, 2839, 10162, 151, 11858, 1579, 2505, 5906, 52, 3174, 1323, 2766, 3336, 6055, 6415, 677, 3445, 7509, 4698, 5057, 12097, 10968, 10240, 4912, 5241, 9369, 3127, 4169, 3482, 787, 6821, 11279, 12231, 241, 11286, 3532, 11404, 6008, 10333, 7280, 2844, 3438, 8077, 975, 5681, 8812, 142, 1105, 4080, 421, 3602, 6221, 4624, 6212, 3263, 8689, 5886, 4782, 5594, 3029, 4213, 504, 605, 9987, 2033, 8291, 10367, 8410, 11316, 11035, 10930, 5435, 3710, 6196, 6950, 5446, 8301, 468, 11973, 11907, 6152, 4948, 11889, 10561, 6153, 6427, 3643, 5415, 56, 9090, 5206, 6760, 1702, 10302, 11635, 3565, 5315, 8214, 7373, 4324, 10120, 11767, 5079, 3262, 11011, 2344, 6715, 1973, 5925, 1018, 3514, 11248, 7500, 7822, 5537, 4749, 8500, 12142, 5456, 7840, 6844, 8429, 7753, 1050, 6118, 3818, 9606, 1190, 5876, 2281, 2031, 5333, 8298, 8320, 12133, 2767, 453, 6381, 418, 3772, 5429, 4774, 1293, 7552, 2361, 1843, 9259, 4115, 218, 2908, 8855, 8760, 2882, 10484, 1954, 2051, 2447, 6147, 576, 3963, 1858, 7535, 3315, 11863, 2925, 347, 3757, 1975, 10596, 3009, 174, 11566, 9551, 5868, 2655, 6554, 1512, 11939, 5383, 10474, 9087, 7796, 6920, 10232, 6374, 1483, 49, 11026, 1489, 2500, 10706, 5942, 1404, 11964, 11143, 948, 4049, 3728, 1159, 5990, 652, 5766, 6190, 11994, 4016, 4077, 2919, 3762, 6328, 7183, 10695, 1962, 7991, 8960, 12121, 9597, 7105, 1200, 6122, 9734, 3956, 1360, 6119, 5297, 3054, 6803, 9166, 1747, 5919, 4433, 3834, 5257, 683, 2459, 8633, 12225, 9786, 9341, 6507, 1566, 11454, 6224, 3570, 8049, 3150, 1319, 4046, 11580, 1958, 7967, 2078, 1112, 11231, 8210, 11367, 441, 1826, 9363, 9118, 4489, 3708, 3238, 11153, 3449, 7080, 1092, 3359, 3205, 8024, 8611, 10361, 11825, 2068, 10900, 4404, 346, 3163, 8257, 7449, 6127, 12164, 11749, 10763, 4222, 8051, 11677, 8921, 8062, 7228, 11071, 11851, 3515, 9011, 5993, 6877, 8080, 1536, 10568, 4103, 9860, 11572, 8700, 1373, 2982, 3448, 11946, 4538, 1908, 4727, 11081, 1866, 7078, 10179, 716, 10125, 6873, 1705, 2450, 11475, 416, 10224, 5826, 7725, 8794, 1756, 4145, 8755, 8328, 5063, 4176, 8524, 10771, 2461, 2275, 8022, 5653, 6693, 6302, 11710, 3889, 212, 6323, 9175, 2769, 5734, 1176, 5508, 11014, 4860, 11164, 11158, 10844, 11841, 1014, 7508, 7365, 10962, 3607, 5232, 8347, 12221, 10029, 7723, 5836, 3200, 1535, 9572, 60, 7784, 10032, 10872, 5676, 3087, 6454, 7406, 3975, 7326, 8545, 2528, 3056, 5845, 5588, 11877, 5102, 1255, 506, 10897, 5784, 9615, 2212, 3338, 9013, 1178, 9513, 6811, 8778, 10347, 3408, 1165, 2575, 10453, 425, 11897, 10104, 377, 4578, 375, 1620, 1038, 11366, 6085, 4167, 6092, 2231, 2800, 12096, 1522, 2151, 8946, 8170, 5002, 12269, 7681, 5163, 10545, 1314, 2894, 3654, 11951, 3947, 9834, 6599, 7350, 7174, 1248, 2442, 8330, 6492, 6330, 10141, 5724, 10964, 1945, 1029, 8945, 6691, 10397, 3624, 6825, 4906, 4670, 512, 7735, 11295, 9389, 12050, 1804, 1403, 6195, 7100, 406, 10602, 7021, 12143, 8914, 9998, 7954, 3393, 8464, 8054, 7376, 8761, 11667, 1737, 4499, 5672, 8307, 9342, 11653, 5609, 4605, 2689, 180, 8151, 5219, 1409, 204, 6780, 9806, 2054, 1344, 9247, 463, 8882, 3981, 1468, 4475, 7043, 3017, 1236, 9168, 4705, 2600, 11232, 4739, 4251, 1226, 6771, 11925, 2360, 3028, 5216, 11839, 10345, 11711, 5368, 11779, 7628, 2622, 6903, 8929, 7605, 7154, 12226, 8481, 8619, 2373, 7302, 10891, 9199, 826, 5043, 5789, 8787, 6671, 10631, 9224, 1506, 7806, 5703, 4719, 11538, 6389, 11379, 4693, 9951, 11872, 9996, 6138, 8820, 4443, 8871, 7186, 10398, 1802, 10734, 1590, 4411, 1223, 2334, 2946, 6828, 2637, 4510, 881, 365, 10362, 1015, 7250, 6742, 2485, 904, 24, 10918, 11009, 11675, 980, 11607, 5082, 7699, 5207, 8239, 844, 7087, 3221, 8016, 8452, 2595, 5289, 6627, 567, 2941, 1406, 2633, 6940, 2945, 3232, 11996, 3769, 7434, 3944, 8190, 6759, 5604, 11024, 9282, 10118, 8809, 9169, 6184, 6643, 6086, 8753, 5370, 8348, 8536, 1282, 3572, 9457, 2021, 4730, 3229, 1706, 3929, 5054, 3154, 9004, 7929, 12282, 1936, 8566, 11444, 11520, 5526, 50, 216, 767, 3805, 4153, 10076, 1279, 11424, 9617, 5170, 12100, 3116, 10080, 1763, 3815, 1734, 1350, 5832, 8420, 4423, 1530, 1694, 10036, 10421, 9559, 5411, 4820, 1160, 9195, 7771, 2840, 9811, 4194, 9270, 7315, 4565, 7211, 10506, 944, 7519, 7002, 8620, 7624, 6883, 3020, 5673, 5410, 1251, 10499, 7014, 2035, 11249, 6164, 10407, 8176, 12217, 10447, 3840, 2712, 4834, 2828, 4352, 1241, 4378, 3451, 4094, 3045, 5781, 9646, 11194, 7592, 8711, 8823, 10588, 7785, 11511, 2626, 530, 10808, 9332, 9349, 2046, 8972, 9757, 8957, 12150, 3268, 3795, 1849, 6513, 4523, 4301, 457, 8, 8835, 3758, 8071, 4390, 10013, 982, 2593, 879, 9687, 10388, 11787, 7171, 6063, 8496, 8443, 1573, 5969, 4649, 9360, 6026, 1030, 11823, 10608, 8468, 11415, 9988, 5650, 12119, 648, 12139, 2307, 8000, 11498, 9855, 9416, 2827, 9754, 11169, 21, 6481];
   uint256[1024] psi_inv_rev = [uint256(1), 10810, 5146, 8246, 11567, 10984, 8155, 6553, 3621, 2545, 3542, 8785, 3195, 5860, 11077, 10643, 9283, 9545, 2975, 11726, 7468, 2639, 2625, 949, 3328, 5777, 7311, 1351, 5828, 7266, 5728, 7698, 4805, 8736, 1062, 2294, 8577, 9154, 7443, 2747, 9970, 1170, 11334, 11499, 3014, 3201, 1326, 5086, 8034, 1177, 2768, 10654, 12149, 10436, 4611, 726, 3051, 9923, 7393, 2963, 9198, 81, 7969, 11289, 8652, 8830, 145, 6747, 8357, 2731, 5911, 7399, 4231, 9741, 8907, 355, 5179, 8595, 8582, 1759, 7901, 1260, 6534, 7657, 9863, 11955, 1428, 1696, 729, 3241, 3289, 2013, 2089, 7197, 9408, 9005, 11462, 6522, 8541, 953, 7222, 10092, 2476, 118, 3949, 8993, 7837, 9893, 12159, 7935, 6915, 9452, 3584, 8112, 9764, 10908, 5331, 4989, 4278, 1673, 480, 2842, 12280, 1022, 9821, 339, 6498, 11745, 10146, 11224, 7644, 404, 11121, 7012, 11082, 3248, 6845, 2381, 4096, 493, 10377, 1378, 4337, 435, 2422, 6250, 2566, 2187, 8643, 9852, 2987, 6267, 8005, 7201, 1002, 5011, 8509, 11414, 1607, 7313, 1067, 7188, 9888, 11847, 3833, 8511, 773, 11899, 4861, 11935, 6591, 9377, 5012, 9808, 9430, 1045, 27, 9223, 3763, 1440, 5084, 1632, 11272, 4885, 11744, 7270, 9611, 3704, 242, 10752, 4143, 4714, 2865, 2370, 8779, 5332, 8311, 9320, 2686, 9042, 2249, 4048, 9405, 1153, 10659, 2126, 5407, 3186, 2399, 3400, 7098, 9153, 671, 3000, 12046, 3016, 10745, 10111, 5559, 11869, 8758, 11813, 4905, 8304, 2645, 8236, 7247, 9984, 7394, 1484, 2780, 5195, 4414, 9442, 4372, 10115, 8232, 3271, 1689, 8925, 113, 4919, 3915, 10123, 4437, 3, 12129, 3149, 3636, 4938, 6998, 9585, 4654, 10863, 10512, 10626, 11848, 922, 4079, 1058, 11177, 10211, 4322, 10331, 709, 8243, 10970, 9139, 4240, 8719, 6065, 835, 10723, 5782, 2948, 2503, 64, 3656, 9830, 11606, 7032, 8455, 7856, 6370, 10542, 3123, 5486, 9235, 6992, 6170, 10929, 8333, 2555, 6167, 11089, 5184, 2692, 168, 3329, 4298, 10327, 1594, 5106, 5961, 8527, 9370, 8212, 8273, 295, 6099, 6523, 11637, 6299, 11130, 8561, 8240, 11341, 1146, 325, 10885, 6347, 1583, 9789, 10800, 1263, 12240, 10806, 5915, 2057, 5369, 4493, 3202, 1815, 6906, 350, 10777, 5735, 9634, 6421, 2738, 723, 12115, 9280, 1693, 10314, 8532, 11942, 9364, 426, 8974, 4754, 10431, 8326, 11713, 6142, 9842, 10238, 10335, 1805, 9407, 3529, 3434, 9381, 12071, 8174, 3030, 10446, 9928, 4737, 10996, 7515, 6860, 8517, 11871, 5908, 11836, 9522, 156, 3969, 3991, 6956, 10258, 10008, 6413, 11099, 2683, 8471, 6171, 11239, 4536, 3860, 5445, 4449, 6833, 147, 3789, 7540, 6752, 4467, 4789, 1041, 8775, 11271, 6364, 10316, 5574, 9945, 1278, 9027, 7210, 522, 2169, 7965, 4916, 4075, 6974, 8724, 654, 1987, 10587, 5529, 7083, 3199, 12233, 6874, 8646, 5862, 6136, 1728, 400, 7341, 6137, 382, 316, 11821, 3988, 6843, 5339, 6093, 8579, 6854, 1359, 1254, 973, 3879, 1922, 3998, 10256, 2302, 11684, 11785, 8076, 9260, 6695, 7507, 6403, 3600, 9026, 6077, 7665, 6068, 8687, 11868, 8209, 11184, 12147, 3477, 6608, 11314, 4212, 8851, 9445, 5009, 1956, 6281, 885, 8757, 1003, 12048, 58, 1010, 5468, 11502, 8807, 8120, 9162, 2920, 7048, 7377, 2049, 1321, 192, 7232, 7591, 4780, 8844, 11612, 5874, 6234, 8953, 9523, 10966, 9115, 12237, 6383, 9784, 10710, 431, 12138, 2127, 9450, 8332, 5808, 12268, 1120, 2535, 9462, 2873, 2434, 791, 4289, 9982, 150, 11641, 170, 6639, 2301, 874, 3821, 1681, 466, 11259, 6263, 2929, 7640, 6320, 10716, 3846, 3793, 6226, 5118, 502, 1901, 2602, 11410, 9696, 11307, 2276, 7899, 4218, 8531, 3454, 12281, 11832, 7988, 7766, 5776, 10440, 8494, 9021, 139, 3332, 2532, 3317, 10243, 2940, 2957, 1481, 11759, 9663, 778, 4504, 1701, 3466, 3578, 4697, 1095, 2643, 6508, 9244, 8195, 8838, 7911, 11048, 7937, 9461, 7455, 9577, 8449, 1842, 72, 4113, 1882, 6125, 1040, 10254, 5275, 1790, 11038, 6879, 6616, 9269, 5406, 4665, 3669, 5287, 4770, 11345, 1783, 5078, 7724, 4974, 3019, 8095, 2478, 9449, 4518, 3094, 11129, 7469, 6878, 2730, 1868, 2253, 10595, 10759, 7866, 3869, 6457, 10939, 10555, 8474, 10526, 2209, 9173, 189, 7119, 2672, 865, 11010, 2213, 8136, 8484, 11522, 12073, 12239, 6763, 769, 845, 3723, 10353, 7, 4360, 3285, 9135, 7235, 8360, 10583, 9060, 7559, 10268, 2832, 8717, 11007, 3753, 3941, 6919, 3536, 6203, 5646, 6105, 3120, 3480, 2171, 3007, 1265, 6685, 5530, 4099, 8345, 4855, 8520, 293, 9057, 9344, 5349, 9656, 10883, 9348, 11722, 5662, 7000, 9694, 3837, 4273, 9068, 5202, 11445, 4050, 7082, 4590, 7207, 682, 11309, 614, 1280, 1371, 12265, 11385, 9804, 5547, 5039, 11274, 1927, 11924, 11408, 7779, 9652, 5461, 9343, 9955, 11066, 7878, 10699, 1555, 10487, 1891, 5103, 3418, 7846, 3469, 6151, 2293, 417, 2338, 7596, 910, 5900, 751, 7570, 6586, 4483, 10783, 3065, 1658, 5618, 3502, 6500, 7246, 11463, 3090, 1398, 4987, 9916, 3670, 3808, 63, 5135, 4684, 3360, 5386, 9667, 4661, 510, 6921, 578, 1944, 450, 7073, 9261, 9929, 364, 5518, 11063, 8038, 7550, 1057, 9689, 7584, 3121, 11053, 9272, 5246, 7814, 10821, 8308, 3407, 11826, 3042, 10945, 10235, 2483, 5509, 12085, 10880, 7070, 4138, 12109, 9600, 7684, 6680, 636, 2947, 3982, 6617, 7790, 10552, 622, 3528, 4913, 4235, 3825, 8896, 4335, 2291, 3375, 146, 5268, 1687, 11883, 5189, 6094, 10886, 10485, 239, 2900, 994, 4554, 11777, 7619, 7383, 5464, 8665, 1892, 5598, 3344, 11260, 10344, 1325, 6565, 2148, 5959, 5797, 3959, 9847, 11041, 5115, 4939, 5690, 2455, 8342, 338, 8635, 9395, 10975, 1744, 7126, 4608, 20, 7287, 4119, 3343, 10138, 10767, 193, 9489, 10058, 6197, 8122, 6204, 923, 11251, 10669, 11914, 7711, 11912, 2185, 392, 11864, 1836, 9714, 11124, 8881, 1942, 3511, 5478, 2776, 11111, 3276, 8951, 10077, 2674, 6505, 1392, 11783, 11034, 7187, 412, 6701, 6444, 9233, 9761, 3744, 4963, 8314, 4883, 5835, 9202, 6613, 1417, 2257, 4505, 12229, 2717, 10754, 9089, 6453, 4566, 2260, 68, 3942, 7057, 8682, 1327, 4924, 4781, 11275, 448, 1445, 1131, 1125, 7429, 1275, 6781, 11113, 6555, 9520, 3114, 5966, 12077, 8400, 579, 5987, 5596, 6636, 4267, 10014, 9828, 1518, 3765, 8113, 7226, 3961, 3534, 8144, 10533, 3495, 4564, 6463, 2065, 11873, 814, 9839, 10584, 5416, 2164, 11573, 2110, 5211, 10423, 1208, 7562, 10381, 7751, 343, 8841, 9307, 10916, 3589, 717, 2429, 8186, 1721, 10753, 4209, 5412, 6296, 3278, 8774, 438, 1218, 5061, 4227, 3368, 612, 4238, 8067, 1526, 540, 125, 6162, 4840, 4032, 9126, 11943, 7885, 1389, 10221, 464, 1928, 3678, 4265, 9084, 8930, 11197, 5209, 8840, 1136, 9051, 8581, 7800, 3171, 2926, 10463];


    uint256 constant q = 12289;
    uint256 nm1modq=12265;//512^-1 mod 12289

//Test vector generated by pythonref
 // q = 12289
 // n = 512
 // prettier-ignore
uint256[]  f_12289_512 = [5072, 8090, 3251, 4761, 259, 10870, 1822, 9635, 5563, 5405, 3608, 11161, 294, 3226, 1677, 3462, 1785, 5110, 12000, 9790, 7452, 906, 2811, 5750, 8453, 5570, 5461, 3876, 9416, 2202, 1440, 5014, 11755, 5161, 10296, 12009, 10171, 762, 9795, 5892, 2884, 6733, 9129, 3340, 1423, 1335, 4953, 7012, 6458, 1816, 10030, 7817, 7748, 729, 11021, 5222, 2349, 8945, 3140, 2398, 2129, 10243, 5994, 1556, 9535, 9028, 4646, 11899, 12083, 1198, 10295, 9255, 9749, 5497, 11392, 4049, 10550, 11244, 8646, 596, 7640, 5203, 6058, 10422, 2769, 2784, 7348, 6243, 1093, 6640, 11396, 6726, 2822, 6098, 3015, 12175, 7507, 9182, 8089, 2086, 6267, 1242, 8082, 2546, 4720, 4335, 10936, 368, 2080, 7503, 11523, 2520, 1058, 8929, 773, 5836, 11288, 696, 5564, 10372, 1651, 10698, 3574, 1320, 5700, 4578, 10861, 2374, 10233, 6116, 318, 5994, 5322, 5558, 5276, 5494, 8935, 5439, 4989, 10624, 1941, 6114, 10181, 5860, 10009, 10573, 919, 3650, 5539, 12058, 3204, 4697, 2973, 6653, 3769, 11981, 6095, 1801, 4735, 281, 2774, 896, 6530, 4436, 2023, 8922, 638, 5135, 9951, 2169, 138, 3770, 2404, 596, 5005, 43, 2811, 9260, 361, 5893, 4499, 12094, 4258, 1530, 10634, 8281, 3924, 11571, 10972, 2128, 8886, 5940, 4747, 6567, 4002, 8715, 11384, 2133, 910, 3965, 12200, 10022, 11833, 2774, 7570, 10516, 2318, 7897, 754, 11516, 7309, 677, 4652, 3835, 4237, 7288, 10272, 11569, 5574, 496, 5151, 2371, 10489, 11654, 5160, 11469, 2002, 5790, 7262, 4514, 4185, 1678, 11862, 1077, 11158, 9803, 7047, 5967, 518, 8277, 6452, 7011, 10680, 10674, 2879, 8136, 2122, 9113, 11702, 12127, 10159, 10528, 6382, 268, 1662, 11568, 10474, 1650, 5405, 8010, 158, 10765, 7669, 5126, 378, 4824, 2246, 6048, 6482, 3401, 3533, 1227, 9056, 10721, 6168, 5591, 1709, 5382, 10322, 119, 10952, 7398, 1089, 8054, 5386, 169, 6123, 10320, 6528, 835, 7537, 9583, 8088, 3305, 2600, 253, 1600, 3917, 5067, 2609, 11206, 9971, 3542, 5237, 8032, 9302, 11116, 8318, 5113, 5532, 4116, 6302, 4835, 10887, 10584, 3121, 10898, 2566, 862, 3682, 4542, 4687, 11756, 7591, 8933, 11449, 7425, 11949, 12144, 12150, 11117, 3893, 6760, 2474, 5608, 4600, 5257, 10302, 9479, 3669, 10474, 6024, 5116, 1694, 10068, 8696, 2952, 7950, 9815, 10984, 1793, 1432, 4982, 2747, 11624, 5523, 12102, 11454, 9363, 11238, 6786, 89, 11561, 5191, 8371, 9858, 2443, 9511, 6935, 3466, 2413, 7184, 2220, 7969, 8114, 2676, 8155, 4700, 1526, 6990, 6560, 2175, 5996, 8774, 7979, 5700, 7978, 9105, 5002, 6389, 3717, 9257, 11191, 5401, 3788, 650, 6862, 2963, 5464, 2379, 9247, 8596, 378, 6691, 10332, 5180, 11576, 2738, 4240, 3689, 10220, 11791, 1200, 589, 106, 1977, 84, 3611, 10622, 4478, 5077, 8266, 8542, 9452, 7064, 195, 2368, 10009, 5252, 4472, 7469, 2963, 2665, 7070, 7021, 11234, 4248, 8827, 3998, 3065, 2318, 8337, 6739, 9458, 10469, 6657, 6753, 1705, 5488, 9437, 9681, 9990, 11040, 2338, 10791, 3418, 3245, 9798, 6672, 1776, 9627, 7953, 3630, 1393, 10813, 10636, 11431, 11199, 4470, 6058, 4495, 8943, 3211, 10085, 10482, 6176, 3616, 1869, 8294, 827, 678, 10921, 2834, 8347, 8421, 4431, 2225, 161, 12040, 6115, 1045, 10236, 9898, 6315, 11291, 2307, 7319, 1994, 7858, 7642, 3512, 3911, 8390, 8500, 6477, 11285, 1243, 9359, 7980, 2372, 2816, 7228];
 // prettier-ignore
uint256[]  g_12289_512 = [2431, 11813, 5051, 9806, 9632, 11571, 5173, 5726, 10072, 11143, 7171, 43, 10749, 8559, 3294, 8255, 10485, 8264, 2616, 10346, 4888, 10482, 1907, 7012, 5565, 12129, 3888, 10224, 1180, 2587, 11320, 11946, 9218, 4278, 1083, 230, 3624, 10027, 8786, 10061, 3740, 5832, 6855, 10326, 3950, 918, 6801, 11512, 9302, 5564, 1717, 1760, 4582, 12271, 2764, 1352, 8818, 4995, 167, 1729, 3579, 11930, 3072, 3546, 6949, 3892, 5740, 4892, 2238, 6937, 2893, 9885, 8881, 5247, 5956, 4764, 1772, 7782, 6748, 177, 6407, 5871, 9708, 8420, 5551, 7744, 8765, 3947, 3512, 4387, 9232, 2560, 3243, 8666, 9734, 10211, 5927, 8535, 4077, 2370, 7915, 10349, 11187, 8490, 7794, 9593, 11754, 5137, 2926, 1144, 11157, 949, 11033, 5779, 1739, 5944, 11952, 10975, 1600, 6043, 8255, 9821, 3059, 4538, 4528, 3521, 6065, 11881, 6248, 10178, 8084, 11550, 6424, 4737, 9073, 9741, 11873, 10752, 10608, 4849, 35, 8181, 9251, 5106, 4213, 6026, 2426, 1704, 3914, 1268, 4908, 4947, 4791, 9632, 2322, 8574, 11258, 5573, 1313, 8879, 12178, 5129, 2487, 4309, 10705, 1660, 9292, 12229, 948, 11717, 3489, 7565, 2902, 5169, 10298, 2582, 4027, 5594, 8625, 5750, 3715, 4013, 11617, 5493, 11283, 4895, 7039, 4159, 12058, 1340, 1847, 10898, 12230, 8001, 3157, 9307, 1522, 10888, 6541, 2413, 5404, 5222, 10867, 6103, 4534, 1913, 6412, 487, 2108, 5039, 7417, 7249, 7443, 4428, 926, 7057, 8092, 1060, 11446, 5202, 10810, 3555, 8742, 9171, 11536, 2875, 10834, 5328, 847, 1896, 314, 11001, 9834, 2864, 4675, 11480, 9327, 2866, 1723, 10899, 5200, 11622, 8694, 2863, 4405, 9019, 9456, 3438, 841, 10431, 7153, 9165, 7181, 473, 6060, 11757, 3790, 879, 10358, 7766, 7563, 1782, 4383, 6526, 9651, 7017, 2276, 10146, 7420, 10175, 9936, 3117, 8828, 583, 8481, 2508, 6454, 7217, 2742, 11316, 9688, 1564, 5503, 7090, 8834, 5455, 876, 342, 7053, 12017, 3137, 801, 11360, 10519, 5985, 10673, 1517, 10327, 664, 5546, 11784, 6397, 8555, 11735, 2135, 2082, 10938, 10131, 2257, 11189, 10687, 4080, 2110, 7692, 9521, 8898, 7391, 7227, 477, 7499, 215, 7208, 5572, 3627, 4325, 2896, 9265, 73, 1944, 10021, 11510, 9753, 2230, 6415, 984, 1415, 4560, 5228, 43, 5652, 6820, 5409, 8429, 2309, 7580, 4981, 5337, 1994, 1775, 9565, 2357, 5304, 4963, 3066, 5461, 2381, 1430, 8590, 370, 9389, 10868, 9710, 356, 10601, 7022, 2242, 3007, 937, 4616, 5150, 1841, 3425, 11406, 2109, 8157, 1745, 11784, 5489, 9564, 4862, 9805, 9696, 6426, 1453, 5988, 7544, 12202, 4649, 9261, 11979, 966, 12088, 11329, 8685, 9553, 5282, 10516, 7427, 11904, 8303, 1040, 7137, 5530, 3374, 5754, 2362, 11468, 11054, 247, 9009, 10892, 11212, 3710, 1675, 9107, 8137, 1060, 6590, 1439, 8653, 6208, 492, 4678, 2214, 8292, 10340, 2461, 2543, 4633, 6684, 5741, 11462, 6209, 6900, 6239, 7672, 9937, 7263, 9113, 7098, 5246, 3612, 1641, 9263, 9407, 2255, 6789, 9329, 1926, 309, 5915, 5001, 473, 2757, 2976, 6151, 11198, 11998, 9986, 2697, 6544, 11328, 8342, 2579, 2413, 8360, 4855, 8519, 11195, 10532, 8160, 5957, 2853, 3745, 4303, 4385, 9069, 3961, 10126, 7391, 10407, 12283, 5288, 8627, 4598, 5389, 6255, 4599, 10662, 10867, 3484, 10357, 6906, 4763, 27, 10064, 5744, 6682, 9988, 5641, 8990, 196, 3652, 453, 5040, 8124, 8722, 11895, 12058, 1126, 8702, 1182];
 // prettier-ignore
uint256[]  f_ntt_12289_512 = [10406, 8420, 7776, 8647, 10107, 8228, 7297, 10926, 3501, 2935, 5129, 11962, 11156, 7609, 5154, 4684, 6430, 1265, 6357, 3401, 6791, 8683, 4150, 6869, 5276, 11240, 10820, 2599, 3575, 6090, 2327, 7449, 10091, 9704, 11656, 4069, 7297, 7382, 9016, 6444, 2730, 2043, 6340, 10818, 8527, 11287, 1102, 1748, 9946, 6134, 2976, 10774, 10666, 5691, 3048, 6671, 8045, 3556, 1089, 8304, 4370, 11489, 8073, 4755, 4853, 6204, 7435, 3393, 24, 4613, 7291, 6435, 2621, 7546, 5170, 3705, 10596, 11556, 2428, 6800, 2986, 4187, 7548, 7817, 6539, 3005, 7208, 10928, 9103, 7881, 10402, 1551, 9587, 2165, 864, 632, 2952, 4156, 10531, 5509, 8327, 9065, 8828, 4986, 9864, 9108, 73, 4106, 11242, 9938, 4492, 8637, 11733, 3703, 2008, 2753, 117, 6927, 2865, 6034, 7921, 8773, 9773, 10340, 8925, 738, 946, 4220, 11973, 5415, 4128, 2641, 1471, 7941, 490, 2136, 6012, 11934, 8917, 2605, 10911, 11443, 928, 6518, 3974, 4893, 3891, 2259, 11676, 6130, 329, 255, 11945, 8037, 2198, 315, 11299, 7401, 8709, 10369, 3488, 1900, 11744, 11590, 6646, 1652, 3312, 6933, 8649, 4685, 9678, 4752, 8316, 10381, 9676, 9141, 4635, 9180, 12064, 10873, 8257, 9911, 2382, 4883, 8228, 5847, 11638, 1147, 3037, 2888, 6620, 9557, 4215, 3682, 1695, 2415, 1, 8500, 7772, 4423, 2128, 10125, 10654, 9608, 8433, 5888, 3570, 10442, 1398, 10647, 10178, 8011, 6003, 10642, 6574, 11955, 388, 6790, 2276, 10788, 3681, 6790, 12053, 5536, 8043, 178, 7422, 3232, 9286, 3057, 2121, 6122, 9891, 5716, 331, 9675, 5674, 6883, 8454, 3976, 11542, 7637, 4716, 7188, 8560, 4665, 3266, 3368, 1710, 7479, 7936, 4930, 9127, 9169, 4726, 6618, 11753, 11431, 10184, 1175, 8050, 10558, 10354, 7200, 2334, 7912, 11452, 8292, 8855, 9157, 10078, 9383, 6904, 2353, 12064, 8558, 5560, 2003, 7173, 12123, 5068, 6525, 4572, 5887, 8825, 11302, 7893, 9685, 3244, 8457, 2686, 10712, 3250, 851, 10945, 11544, 4826, 11341, 7494, 2428, 4188, 8754, 2972, 9687, 9913, 676, 156, 10861, 12139, 3158, 165, 2129, 2657, 2755, 2708, 11352, 5258, 7384, 5308, 5440, 10059, 8614, 7221, 10502, 8790, 6124, 3431, 952, 7358, 1631, 3066, 9394, 5138, 7157, 5884, 2687, 10259, 5053, 9079, 6092, 553, 7807, 3130, 5141, 1200, 10240, 6394, 6360, 8349, 9027, 9947, 6386, 615, 2091, 3541, 4074, 3468, 992, 4103, 10906, 5105, 8278, 701, 4334, 1823, 6151, 4983, 4684, 8924, 2912, 10295, 3715, 6561, 3702, 1992, 9406, 2175, 2279, 8603, 7267, 3744, 1442, 534, 4119, 7014, 8352, 7391, 5106, 10284, 9061, 4553, 2421, 10066, 6373, 4668, 5490, 1189, 10698, 2365, 3736, 8108, 5851, 11733, 10118, 694, 7199, 9614, 6897, 1502, 1734, 5064, 11541, 3216, 9593, 7139, 9974, 9886, 3964, 1288, 6216, 8357, 9460, 11530, 898, 4678, 3334, 2415, 9454, 1628, 352, 1588, 9287, 5350, 9565, 104, 9041, 8817, 4233, 9566, 5327, 1208, 4537, 1141, 4929, 5687, 1424, 942, 5266, 5443, 2389, 7962, 469, 11823, 9818, 4662, 8858, 4439, 7273, 9993, 6357, 1814, 7332, 1307, 7394, 10094, 4717, 10879, 3500, 63, 3945, 11660, 8270, 9140, 5362, 12014, 6845, 7187, 9343, 4754, 8062, 8744, 10365, 6450, 4492, 7748, 10390, 9713, 11143, 4370, 5462, 10072, 7694, 1672, 11895, 1263, 8632, 10679, 8134, 9211, 10825, 9325, 2776, 5248, 130, 1364, 3280, 5997, 8689, 1304, 840, 7365, 7458];
 // prettier-ignore
uint256[]  g_ntt_12289_512 = [5530, 6817, 6066, 8405, 147, 4396, 10893, 6991, 8484, 4928, 7318, 3032, 3603, 672, 5414, 4656, 3103, 9045, 10494, 830, 5780, 5879, 7150, 6507, 6874, 1497, 6357, 5084, 7463, 1510, 7934, 91, 4643, 7793, 4687, 5529, 7614, 3923, 3262, 6007, 4064, 9736, 2152, 2522, 3225, 6688, 3822, 4702, 10402, 3076, 6065, 4895, 7290, 2605, 3769, 9145, 2961, 770, 5622, 9064, 9410, 3552, 10102, 564, 7696, 11832, 1260, 7290, 1818, 11890, 223, 11055, 2719, 5351, 9194, 4459, 3558, 959, 10238, 6743, 6450, 5180, 4703, 8552, 2152, 10771, 4952, 2542, 12089, 3356, 8873, 12079, 7179, 9629, 7897, 8024, 7246, 1473, 11880, 6340, 2715, 5704, 355, 5859, 3138, 6189, 10097, 2973, 10948, 6056, 11874, 7106, 1318, 10070, 5870, 4922, 6383, 3434, 8125, 2394, 7722, 8255, 9342, 5938, 1392, 6837, 10464, 5497, 2342, 7006, 2852, 4496, 11058, 6491, 10921, 1590, 9039, 2052, 4066, 895, 10939, 9768, 7045, 11711, 1465, 8411, 7339, 9428, 4312, 6811, 4622, 4087, 8713, 9766, 732, 6848, 6600, 2746, 11434, 2502, 3828, 2413, 11434, 4942, 8778, 12153, 9598, 5892, 4599, 1620, 11052, 8423, 3740, 7004, 2261, 11975, 1516, 11924, 2114, 6812, 3287, 52, 5096, 9315, 3731, 5403, 4829, 4906, 1395, 5933, 8453, 25, 6911, 5365, 1132, 4409, 3484, 6557, 1646, 2155, 4739, 7996, 10230, 843, 8451, 1697, 1692, 4719, 10749, 12279, 5773, 8555, 2181, 2907, 11521, 8769, 7405, 4232, 8773, 949, 11856, 6132, 7987, 2343, 8605, 8288, 9128, 2221, 1404, 2501, 5727, 7050, 9024, 7856, 3455, 5846, 3035, 7900, 6300, 10543, 6753, 11516, 2595, 7860, 6182, 83, 12105, 9388, 12151, 12037, 6246, 4636, 6274, 4397, 5476, 5303, 4875, 6049, 1756, 9939, 3752, 10680, 3993, 816, 12127, 5239, 11058, 7451, 11372, 1107, 7702, 6058, 10355, 5116, 6320, 7195, 10611, 3311, 134, 1076, 6881, 5787, 931, 9321, 876, 8318, 5862, 9295, 7651, 10654, 11679, 7384, 9927, 9296, 9874, 3053, 5493, 602, 6805, 8867, 2192, 2194, 3596, 3586, 11038, 8761, 9858, 10537, 6769, 7962, 1887, 8718, 7211, 6636, 8310, 9869, 6601, 9080, 4489, 10686, 9339, 8280, 6903, 6950, 7924, 2011, 8330, 9496, 10373, 11821, 9656, 8963, 7311, 827, 6626, 10688, 11547, 9763, 5189, 2257, 2949, 3726, 7866, 4970, 7356, 4658, 4834, 5505, 8693, 3595, 1768, 4581, 9794, 5243, 2657, 11335, 3007, 474, 2607, 11159, 10553, 5637, 10290, 9916, 4384, 7706, 3560, 3331, 11297, 1151, 4076, 8752, 8065, 9823, 7843, 8940, 1583, 9105, 2094, 10753, 6326, 922, 2393, 4519, 3053, 10874, 6788, 3892, 11501, 9764, 8351, 1001, 316, 11136, 5379, 10174, 2239, 3258, 7824, 2435, 7617, 6715, 9101, 751, 2801, 1864, 9483, 8716, 7002, 8484, 11340, 4815, 1615, 2799, 656, 5059, 10698, 1280, 554, 6293, 4212, 4545, 9811, 11971, 1059, 11821, 4176, 8643, 6426, 8601, 9257, 8915, 3564, 5453, 4940, 11815, 8890, 11710, 1167, 2394, 7327, 3173, 5141, 8803, 1616, 5528, 213, 792, 2889, 4275, 11148, 5469, 4791, 10441, 10593, 6469, 3376, 4891, 3351, 539, 4770, 9302, 11843, 5654, 8679, 8553, 10036, 6854, 200, 5119, 572, 4414, 6398, 8487, 8665, 5235, 4027, 2718, 3932, 8590, 2413, 7300, 2714, 4323, 1977, 5539, 195, 3465, 1924, 10158, 2356, 12029, 3498, 10777, 3482, 11580, 6593, 11284, 6543, 1176, 6602, 9560, 2618, 2571, 355, 1239, 1976, 1444, 9066, 11538, 3318, 12286];
 // prettier-ignore
uint256[]  f_ntt_mul_g_ntt_12289_512 = [8082, 9510, 4034, 889, 11049, 3761, 969, 7531, 12260, 11816, 3416, 3945, 10038, 1024, 7726, 8018, 7243, 866, 5666, 8649, 914, 11140, 6854, 1490, 2385, 2639, 1207, 2641, 806, 3728, 4340, 1964, 6845, 9055, 7067, 8631, 789, 6702, 2615, 11047, 10042, 7046, 2890, 1416, 9082, 8418, 9006, 10044, 9490, 4569, 9188, 6631, 2637, 4521, 9986, 3699, 5163, 9962, 2436, 9620, 2706, 9448, 3642, 2818, 2417, 3531, 3882, 9502, 6765, 2763, 3745, 10193, 11168, 9281, 11417, 4179, 10205, 9815, 9506, 2141, 2837, 10864, 7612, 11113, 1023, 9918, 6760, 5836, 10461, 2708, 6556, 6093, 6673, 4641, 2613, 8100, 7332, 1866, 6260, 1722, 8334, 6937, 245, 2021, 9530, 12058, 12030, 4161, 3081, 5295, 3748, 3256, 4532, 4384, 1809, 7788, 9471, 8103, 2759, 5821, 3609, 2038, 4385, 3076, 11710, 7216, 6299, 7997, 9557, 1347, 194, 2762, 7971, 4965, 5575, 4476, 510, 8880, 3972, 8854, 4661, 6769, 12, 5319, 9213, 11451, 8702, 1015, 11168, 5697, 9091, 9909, 1244, 11788, 11366, 6545, 3748, 9429, 939, 1159, 6210, 903, 11282, 11040, 2705, 8819, 9222, 600, 9547, 7387, 10089, 823, 10670, 6800, 3016, 5352, 9641, 4197, 3621, 1073, 6647, 11523, 9429, 3556, 746, 8611, 2305, 11109, 9199, 3638, 7043, 5434, 4935, 5507, 1656, 5461, 3484, 3885, 12152, 7590, 7612, 11857, 11568, 1093, 3372, 979, 6541, 9197, 9944, 4131, 3885, 10641, 4758, 4881, 1947, 8225, 9803, 3598, 10012, 1075, 3697, 1148, 7574, 5953, 10656, 584, 11048, 1496, 11204, 1799, 5435, 1132, 1377, 890, 728, 6072, 3701, 9164, 11963, 1189, 6288, 7608, 10465, 5147, 1486, 6236, 1217, 11476, 9800, 7798, 6719, 10229, 8347, 8173, 11231, 10159, 4557, 8205, 2609, 3775, 9527, 7865, 3326, 1058, 2851, 171, 10360, 6889, 2994, 10663, 3432, 5589, 5807, 7017, 3524, 6920, 9960, 8162, 2640, 5719, 9015, 8367, 4538, 2342, 919, 11475, 681, 5150, 8353, 10219, 8266, 5404, 4125, 9069, 1464, 11269, 1845, 6887, 9609, 11037, 213, 10858, 8171, 8868, 10727, 11427, 1723, 7189, 4637, 702, 4130, 4232, 1076, 8437, 2321, 6364, 3922, 10225, 11530, 4870, 3885, 10853, 2379, 4529, 10197, 1786, 8305, 7777, 9844, 10899, 1095, 6483, 8734, 7830, 6676, 11552, 7002, 4393, 7194, 10542, 8649, 819, 5713, 1939, 3698, 4311, 1761, 439, 11312, 9105, 737, 6446, 1700, 1325, 7352, 9017, 7204, 3226, 5091, 2087, 10378, 1753, 11936, 1311, 4182, 933, 6453, 7663, 7761, 9104, 7774, 9275, 10320, 1595, 3937, 8302, 2105, 6463, 11297, 8589, 3641, 2312, 12095, 8215, 6304, 3938, 6410, 1239, 6948, 3093, 12226, 2488, 10294, 753, 2745, 1755, 7747, 2480, 8815, 3300, 6411, 1532, 2912, 4016, 2232, 11637, 9760, 8753, 9909, 1323, 11552, 11346, 7882, 11631, 1075, 12121, 1294, 10852, 790, 1401, 3988, 8778, 585, 9372, 1535, 391, 8060, 1361, 3589, 4458, 2472, 2612, 7161, 3429, 9911, 3427, 3888, 6893, 5110, 9145, 2936, 5482, 4028, 9817, 10309, 6912, 4022, 4701, 7196, 816, 9218, 8849, 3992, 7189, 7364, 11084, 5773, 7877, 11307, 10081, 1324, 10603, 6950, 10787, 9834, 12003, 6168, 872, 311, 3628, 8882, 5450, 6658, 1127, 1191, 11140, 1454, 5200, 1159, 4065, 11348, 1127, 5764, 2296, 5702, 823, 1529, 10746, 2204, 10450, 11862, 2667, 11381, 5856, 10593, 12123, 2966, 9804, 2317, 11085, 8049, 6609, 162, 2427, 4949, 8550, 3476, 12136, 46, 8188, 6538, 2204];
 // prettier-ignore
uint256[]  f_mul_g_12289_512 = [7776, 3771, 1701, 793, 6632, 11882, 8190, 10514, 3180, 7201, 330, 8803, 9526, 9628, 11277, 10452, 1154, 8464, 5665, 5922, 1637, 4581, 1463, 8241, 12011, 10363, 3163, 10545, 3246, 926, 7428, 9116, 8447, 1640, 5444, 11091, 8091, 3629, 9925, 577, 3645, 941, 5525, 8925, 2269, 3850, 1236, 3790, 4098, 11075, 4234, 1052, 11576, 11490, 3178, 10762, 9572, 9936, 5496, 71, 2997, 2607, 11191, 4730, 5641, 181, 11524, 7832, 3221, 11617, 252, 7365, 1727, 9798, 11251, 9945, 3372, 9672, 10069, 3010, 12077, 6923, 11497, 10732, 7213, 7937, 8934, 11292, 447, 3018, 11427, 1828, 4012, 10361, 7585, 8663, 2196, 5737, 843, 233, 5916, 7684, 916, 6022, 4379, 6543, 10388, 1399, 8198, 2131, 7159, 3916, 5547, 3826, 934, 9383, 1209, 9993, 8019, 8265, 11953, 11250, 1606, 3959, 2076, 4243, 12113, 1385, 4274, 8752, 167, 6008, 748, 9073, 10282, 426, 9999, 7, 3815, 8096, 11296, 3553, 2728, 8372, 4944, 1728, 4112, 52, 731, 7071, 2704, 3857, 7167, 8642, 8448, 5081, 9581, 5424, 427, 64, 7386, 2991, 7295, 6506, 11520, 4844, 7363, 5023, 2927, 1469, 3772, 2257, 7699, 6933, 9591, 6138, 3886, 10769, 106, 8474, 2726, 11196, 4269, 7760, 6376, 10953, 7725, 3395, 7012, 9702, 1859, 9634, 5960, 10206, 4340, 1314, 6659, 9976, 2715, 10972, 1313, 5979, 9250, 7, 9745, 6230, 9229, 4183, 245, 6585, 9160, 2828, 204, 8349, 12166, 9996, 8112, 11968, 9293, 4768, 10885, 539, 4291, 5781, 8691, 706, 8366, 8499, 3200, 1558, 2774, 9443, 643, 1842, 8263, 8292, 3309, 9343, 3186, 12057, 10924, 3566, 1804, 3495, 7863, 1265, 212, 2594, 11590, 10983, 11807, 6688, 5874, 9910, 9334, 11254, 9853, 7724, 11034, 8143, 130, 5403, 10071, 10058, 7611, 6424, 7397, 2264, 9832, 1678, 11284, 9912, 10957, 4846, 1222, 11189, 2533, 7643, 905, 1378, 10769, 2076, 7250, 11646, 2952, 5248, 3447, 1719, 6660, 7780, 11741, 3295, 5895, 1082, 4722, 7010, 1857, 12126, 2975, 228, 6721, 10300, 932, 9478, 10529, 1245, 3831, 6059, 9701, 1590, 3507, 3305, 3451, 10726, 9501, 3487, 8199, 403, 8105, 2407, 8548, 1400, 10759, 9964, 9505, 2485, 9528, 2907, 10018, 9016, 3756, 5482, 6053, 11380, 4801, 5335, 8221, 2831, 11397, 10583, 6977, 4936, 9673, 8916, 9898, 5051, 7887, 1245, 8586, 6243, 6348, 8975, 6392, 6669, 8362, 1324, 9401, 4164, 8545, 11435, 737, 8424, 4193, 154, 3734, 11801, 10241, 5269, 2123, 3162, 9701, 5213, 1722, 6921, 4121, 3165, 4945, 8964, 6177, 647, 3883, 11273, 9614, 6664, 6704, 10607, 4063, 5115, 4005, 6531, 2542, 373, 2737, 11732, 11303, 8071, 10366, 1923, 3637, 6702, 9366, 2817, 12140, 9156, 8684, 4386, 3790, 5693, 11033, 8785, 190, 1324, 4701, 7922, 4058, 10320, 2394, 8730, 2279, 10712, 9852, 8781, 122, 91, 10340, 8507, 3391, 7056, 3236, 350, 11294, 7428, 5790, 2474, 121, 11276, 10097, 133, 2242, 5244, 7153, 6836, 10487, 3140, 3507, 6679, 2056, 9955, 1293, 264, 512, 9551, 3768, 3578, 8305, 2291, 7814, 937, 2811, 7189, 3593, 3888, 2592, 10261, 4808, 2086, 4607, 3874, 6660, 8353, 5204, 12015, 8950, 8541, 6971, 12114, 1670, 1998, 7940, 9703, 9283, 2477, 10308, 12165, 6611, 4638, 10324, 953, 2687, 9684, 1981, 8440, 10168, 6522, 9165, 8864, 6198, 3778, 9598, 2403, 9295, 9258, 356, 11895, 2172, 10951, 4184, 3503, 1777, 2002, 10839, 8829];

/*
def ntt(self, f):
        # following eprint 2016/504 Algorithm 1
        a = [_ for _ in f]
        n = len(a)
        t = n
        m = 1
        while m < n:
            t //= 2
            for i in range(m):
                j1 = 2*i*t
                j2 = j1+t-1
                S = self.ψ_rev[m+i]
                for j in range(j1, j2+1):
                    U = a[j]
                    V = a[j+t]*S
                    a[j] = (U+V) % self.q
                    a[j+t] = (U-V) % self.q//end for j, end for i
            m = 2*m
        return a
*/


    function test_NWCntt() public view{
        console.log("length:",f_12289_512.length);
        uint256[] memory g=ntt.NWC_ntt(f_12289_512);
        for(uint i=0;i<512;i++)
          assertEq(f_ntt_12289_512[i], g[i]);
    }


    function testbench_NWCntt() public view{
        uint256[] memory g=ntt.NWC_ntt(f_12289_512);
    }

    //test inverse NTT
    function test_NWCIntt() public view{
        console.log("length:",f_ntt_mul_g_ntt_12289_512.length);
        uint256[] memory g=ntt.NWC_Intt(f_ntt_mul_g_ntt_12289_512);
        for(uint i=0;i<512;i++)
           assertEq(f_mul_g_12289_512[i], g[i]);
    }


    function testbench_NWCIntt() public view{
        uint256[] memory g=ntt.NWC_Intt(f_12289_512);
    }


    function test_mulNTTPoly() public view{
        uint256[] memory g=ntt.mul_NTTPoly(f_12289_512, g_12289_512);
        for(uint i=0;i<512;i++)
           assertEq(f_mul_g_12289_512[i], g[i]);
        }

    
    function test_mul_halfNTTPoly() public view{
        uint256[] memory g=ntt.mul_halfNTTPoly(f_12289_512, g_ntt_12289_512);
        for(uint i=0;i<512;i++)
           assertEq(f_mul_g_12289_512[i], g[i]);
        }

    function testbench_halfNTTPoly() public view{
        uint256[] memory g=ntt.mul_halfNTTPoly(f_12289_512, g_ntt_12289_512);
        } 

}