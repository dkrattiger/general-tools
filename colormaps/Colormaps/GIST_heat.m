function cm_data=GIST_heat(m)

cm =[...
     0     0     0
     1     0     0
     1     0     0
     2     0     0
     2     0     0
     4     0     0
     4     0     0
     5     0     0
     5     0     0
     7     0     0
     7     0     0
     8     0     0
     8     0     0
    10     0     0
    10     0     0
    11     0     0
    11     0     0
    13     0     0
    13     0     0
    15     0     0
    15     0     0
    17     0     0
    17     0     0
    18     0     0
    18     0     0
    20     0     0
    20     0     0
    21     0     0
    21     0     0
    23     0     0
    23     0     0
    24     0     0
    24     0     0
    26     0     0
    26     0     0
    27     0     0
    27     0     0
    28     0     0
    28     0     0
    30     0     0
    30     0     0
    31     0     0
    31     0     0
    33     0     0
    33     0     0
    34     0     0
    34     0     0
    36     0     0
    36     0     0
    37     0     0
    37     0     0
    39     0     0
    39     0     0
    40     0     0
    40     0     0
    42     0     0
    42     0     0
    43     0     0
    43     0     0
    46     0     0
    46     0     0
    47     0     0
    47     0     0
    49     0     0
    49     0     0
    50     0     0
    50     0     0
    52     0     0
    52     0     0
    53     0     0
    53     0     0
    55     0     0
    55     0     0
    56     0     0
    56     0     0
    57     0     0
    57     0     0
    59     0     0
    59     0     0
    60     0     0
    60     0     0
    62     0     0
    62     0     0
    63     0     0
    63     0     0
    65     0     0
    65     0     0
    66     0     0
    66     0     0
    68     0     0
    68     0     0
    69     0     0
    69     0     0
    70     0     0
    70     0     0
    72     0     0
    72     0     0
    73     0     0
    73     0     0
    76     0     0
    76     0     0
    78     0     0
    78     0     0
    79     0     0
    79     0     0
    81     0     0
    81     0     0
    82     0     0
    82     0     0
    84     0     0
    84     0     0
    85     0     0
    85     0     0
    86     0     0
    86     0     0
    88     0     0
    88     0     0
    89     0     0
    89     0     0
    92     0     0
    92     0     0
    94     0     0
    94     0     0
    95     0     0
    95     0     0
    97     0     0
    97     0     0
    98     0     0
    98     0     0
    99     0     0
    99     0     0
   101     0     0
   101     0     0
   102     0     0
   102     0     0
   104     0     0
   104     0     0
   105     0     0
   105     0     0
   108     0     0
   108     0     0
   110     0     0
   110     0     0
   111     0     0
   111     0     0
   113     0     0
   113     0     0
   114     0     0
   114     0     0
   115     0     0
   115     0     0
   117     0     0
   117     0     0
   118     0     0
   118     0     0
   120     0     0
   120     0     0
   121     0     0
   121     0     0
   123     0     0
   123     0     0
   124     0     0
   124     0     0
   126     0     0
   126     0     0
   127     0     0
   127     0     0
   128     0     0
   128     0     0
   130     0     0
   130     0     0
   131     0     0
   131     0     0
   133     0     0
   133     0     0
   134     0     0
   134     0     0
   136     0     0
   136     0     0
   139     0     0
   139     0     0
   140     0     0
   140     0     0
   141     0     0
   141     0     0
   143     0     0
   143     0     0
   144     0     0
   144     0     0
   146     0     0
   146     0     0
   147     0     0
   147     0     0
   149     0     0
   149     0     0
   150     0     0
   150     0     0
   152     0     0
   152     0     0
   153     0     0
   153     0     0
   155     0     0
   155     0     0
   156     0     0
   156     0     0
   157     0     0
   157     0     0
   159     0     0
   159     0     0
   160     0     0
   160     0     0
   162     0     0
   162     0     0
   163     0     0
   163     0     0
   165     0     0
   165     0     0
   166     0     0
   166     0     0
   169     0     0
   169     0     0
   170     0     0
   170     0     0
   172     0     0
   172     0     0
   173     0     0
   173     0     0
   175     1     0
   175     1     0
   176     3     0
   176     3     0
   178     5     0
   178     5     0
   179     7     0
   179     7     0
   181     9     0
   181     9     0
   182    11     0
   182    11     0
   185    15     0
   185    15     0
   186    17     0
   186    17     0
   188    18     0
   188    18     0
   189    20     0
   189    20     0
   191    22     0
   191    22     0
   192    24     0
   192    24     0
   194    26     0
   194    26     0
   195    28     0
   195    28     0
   197    30     0
   197    30     0
   198    32     0
   198    32     0
   201    35     0
   201    35     0
   202    37     0
   202    37     0
   204    39     0
   204    39     0
   205    41     0
   205    41     0
   207    43     0
   207    43     0
   208    45     0
   208    45     0
   210    47     0
   210    47     0
   211    49     0
   211    49     0
   212    51     0
   212    51     0
   214    52     0
   214    52     0
   215    54     0
   215    54     0
   217    56     0
   217    56     0
   218    58     0
   218    58     0
   220    60     0
   220    60     0
   221    62     0
   221    62     0
   223    64     0
   223    64     0
   224    66     0
   224    66     0
   226    68     0
   226    68     0
   227    69     0
   227    69     0
   228    71     0
   228    71     0
   231    75     0
   231    75     0
   233    77     0
   233    77     0
   234    79     0
   234    79     0
   236    81     0
   236    81     0
   237    83     0
   237    83     0
   239    85     0
   239    85     0
   240    86     0
   240    86     0
   241    88     0
   241    88     0
   243    90     0
   243    90     0
   244    92     0
   244    92     0
   246    94     0
   246    94     0
   247    96     0
   247    96     0
   249    98     0
   249    98     0
   250   100     0
   250   100     0
   252   102     0
   252   102     0
   253   103     0
   253   103     0
   255   105     0
   255   105     0
   255   107     0
   255   107     0
   255   109     0
   255   109     0
   255   111     0
   255   111     0
   255   115     0
   255   115     0
   255   117     0
   255   117     0
   255   119     0
   255   119     0
   255   120     0
   255   120     0
   255   122     0
   255   122     0
   255   124     0
   255   124     0
   255   126     0
   255   126     0
   255   128     0
   255   128     0
   255   130     0
   255   130     0
   255   132     0
   255   132     0
   255   136     7
   255   136     7
   255   137    11
   255   137    11
   255   139    15
   255   139    15
   255   141    19
   255   141    19
   255   143    23
   255   143    23
   255   145    27
   255   145    27
   255   147    31
   255   147    31
   255   149    35
   255   149    35
   255   151    39
   255   151    39
   255   153    43
   255   153    43
   255   156    51
   255   156    51
   255   158    54
   255   158    54
   255   160    58
   255   160    58
   255   162    62
   255   162    62
   255   164    66
   255   164    66
   255   166    70
   255   166    70
   255   168    74
   255   168    74
   255   170    78
   255   170    78
   255   171    82
   255   171    82
   255   173    86
   255   173    86
   255   175    90
   255   175    90
   255   177    94
   255   177    94
   255   179    98
   255   179    98
   255   181   102
   255   181   102
   255   183   105
   255   183   105
   255   185   109
   255   185   109
   255   187   113
   255   187   113
   255   188   117
   255   188   117
   255   190   121
   255   190   121
   255   192   125
   255   192   125
   255   196   133
   255   196   133
   255   198   137
   255   198   137
   255   200   141
   255   200   141
   255   202   145
   255   202   145
   255   204   149
   255   204   149
   255   205   153
   255   205   153
   255   207   156
   255   207   156
   255   209   160
   255   209   160
   255   211   164
   255   211   164
   255   213   168
   255   213   168
   255   215   172
   255   215   172
   255   217   176
   255   217   176
   255   219   180
   255   219   180
   255   221   184
   255   221   184
   255   222   188
   255   222   188
   255   224   192
   255   224   192
   255   226   196
   255   226   196
   255   228   200
   255   228   200
   255   230   204
   255   230   204
   255   232   207
   255   232   207
   255   236   215
   255   236   215
   255   238   219
   255   238   219
   255   239   223
   255   239   223
   255   241   227
   255   241   227
   255   243   231
   255   243   231
   255   245   235
   255   245   235
   255   247   239
   255   247   239
   255   249   243
   255   249   243
   255   251   247
   255   251   247
   255   253   251]/255;
% A = [0.000000e+00   0   0   0 1.000000e+00   1   0   0
% 1.000000e+00   1   0   0 2.000000e+00   2   0   0
% 2.000000e+00   2   0   0 3.000000e+00   4   0   0
% 3.000000e+00   4   0   0 4.000000e+00   5   0   0
% 4.000000e+00   5   0   0 5.000000e+00   7   0   0
% 5.000000e+00   7   0   0 6.000000e+00   8   0   0
% 6.000000e+00   8   0   0 7.000000e+00  10   0   0
% 7.000000e+00  10   0   0 8.000000e+00  11   0   0
% 8.000000e+00  11   0   0 9.000000e+00  13   0   0
% 9.000000e+00  13   0   0 1.000000e+01  15   0   0
% 1.000000e+01  15   0   0 1.100000e+01  17   0   0
% 1.100000e+01  17   0   0 1.200000e+01  18   0   0
% 1.200000e+01  18   0   0 1.300000e+01  20   0   0
% 1.300000e+01  20   0   0 1.400000e+01  21   0   0
% 1.400000e+01  21   0   0 1.500000e+01  23   0   0
% 1.500000e+01  23   0   0 1.600000e+01  24   0   0
% 1.600000e+01  24   0   0 1.700000e+01  26   0   0
% 1.700000e+01  26   0   0 1.800000e+01  27   0   0
% 1.800000e+01  27   0   0 1.900000e+01  28   0   0
% 1.900000e+01  28   0   0 2.000000e+01  30   0   0
% 2.000000e+01  30   0   0 2.100000e+01  31   0   0
% 2.100000e+01  31   0   0 2.200000e+01  33   0   0
% 2.200000e+01  33   0   0 2.300000e+01  34   0   0
% 2.300000e+01  34   0   0 2.400000e+01  36   0   0
% 2.400000e+01  36   0   0 2.500000e+01  37   0   0
% 2.500000e+01  37   0   0 2.600000e+01  39   0   0
% 2.600000e+01  39   0   0 2.700000e+01  40   0   0
% 2.700000e+01  40   0   0 2.800000e+01  42   0   0
% 2.800000e+01  42   0   0 2.900000e+01  43   0   0
% 2.900000e+01  43   0   0 3.000000e+01  46   0   0
% 3.000000e+01  46   0   0 3.100000e+01  47   0   0
% 3.100000e+01  47   0   0 3.200000e+01  49   0   0
% 3.200000e+01  49   0   0 3.300000e+01  50   0   0
% 3.300000e+01  50   0   0 3.400000e+01  52   0   0
% 3.400000e+01  52   0   0 3.500000e+01  53   0   0
% 3.500000e+01  53   0   0 3.600000e+01  55   0   0
% 3.600000e+01  55   0   0 3.700000e+01  56   0   0
% 3.700000e+01  56   0   0 3.800000e+01  57   0   0
% 3.800000e+01  57   0   0 3.900000e+01  59   0   0
% 3.900000e+01  59   0   0 4.000000e+01  60   0   0
% 4.000000e+01  60   0   0 4.100000e+01  62   0   0
% 4.100000e+01  62   0   0 4.200000e+01  63   0   0
% 4.200000e+01  63   0   0 4.300000e+01  65   0   0
% 4.300000e+01  65   0   0 4.400000e+01  66   0   0
% 4.400000e+01  66   0   0 4.500000e+01  68   0   0
% 4.500000e+01  68   0   0 4.600000e+01  69   0   0
% 4.600000e+01  69   0   0 4.700000e+01  70   0   0
% 4.700000e+01  70   0   0 4.800000e+01  72   0   0
% 4.800000e+01  72   0   0 4.900000e+01  73   0   0
% 4.900000e+01  73   0   0 5.000000e+01  76   0   0
% 5.000000e+01  76   0   0 5.100000e+01  78   0   0
% 5.100000e+01  78   0   0 5.200000e+01  79   0   0
% 5.200000e+01  79   0   0 5.300000e+01  81   0   0
% 5.300000e+01  81   0   0 5.400000e+01  82   0   0
% 5.400000e+01  82   0   0 5.500000e+01  84   0   0
% 5.500000e+01  84   0   0 5.600000e+01  85   0   0
% 5.600000e+01  85   0   0 5.700000e+01  86   0   0
% 5.700000e+01  86   0   0 5.800000e+01  88   0   0
% 5.800000e+01  88   0   0 5.900000e+01  89   0   0
% 5.900000e+01  89   0   0 6.000000e+01  92   0   0
% 6.000000e+01  92   0   0 6.100000e+01  94   0   0
% 6.100000e+01  94   0   0 6.200000e+01  95   0   0
% 6.200000e+01  95   0   0 6.300000e+01  97   0   0
% 6.300000e+01  97   0   0 6.400000e+01  98   0   0
% 6.400000e+01  98   0   0 6.500000e+01  99   0   0
% 6.500000e+01  99   0   0 6.600000e+01 101   0   0
% 6.600000e+01 101   0   0 6.700000e+01 102   0   0
% 6.700000e+01 102   0   0 6.800000e+01 104   0   0
% 6.800000e+01 104   0   0 6.900000e+01 105   0   0
% 6.900000e+01 105   0   0 7.000000e+01 108   0   0
% 7.000000e+01 108   0   0 7.100000e+01 110   0   0
% 7.100000e+01 110   0   0 7.200000e+01 111   0   0
% 7.200000e+01 111   0   0 7.300000e+01 113   0   0
% 7.300000e+01 113   0   0 7.400000e+01 114   0   0
% 7.400000e+01 114   0   0 7.500000e+01 115   0   0
% 7.500000e+01 115   0   0 7.600000e+01 117   0   0
% 7.600000e+01 117   0   0 7.700000e+01 118   0   0
% 7.700000e+01 118   0   0 7.800000e+01 120   0   0
% 7.800000e+01 120   0   0 7.900000e+01 121   0   0
% 7.900000e+01 121   0   0 8.000000e+01 123   0   0
% 8.000000e+01 123   0   0 8.100000e+01 124   0   0
% 8.100000e+01 124   0   0 8.200000e+01 126   0   0
% 8.200000e+01 126   0   0 8.300000e+01 127   0   0
% 8.300000e+01 127   0   0 8.400000e+01 128   0   0
% 8.400000e+01 128   0   0 8.500000e+01 130   0   0
% 8.500000e+01 130   0   0 8.600000e+01 131   0   0
% 8.600000e+01 131   0   0 8.700000e+01 133   0   0
% 8.700000e+01 133   0   0 8.800000e+01 134   0   0
% 8.800000e+01 134   0   0 8.900000e+01 136   0   0
% 8.900000e+01 136   0   0 9.000000e+01 139   0   0
% 9.000000e+01 139   0   0 9.100000e+01 140   0   0
% 9.100000e+01 140   0   0 9.200000e+01 141   0   0
% 9.200000e+01 141   0   0 9.300000e+01 143   0   0
% 9.300000e+01 143   0   0 9.400000e+01 144   0   0
% 9.400000e+01 144   0   0 9.500000e+01 146   0   0
% 9.500000e+01 146   0   0 9.600000e+01 147   0   0
% 9.600000e+01 147   0   0 9.700000e+01 149   0   0
% 9.700000e+01 149   0   0 9.800000e+01 150   0   0
% 9.800000e+01 150   0   0 9.900000e+01 152   0   0
% 9.900000e+01 152   0   0 1.000000e+02 153   0   0
% 1.000000e+02 153   0   0 1.010000e+02 155   0   0
% 1.010000e+02 155   0   0 1.020000e+02 156   0   0
% 1.020000e+02 156   0   0 1.030000e+02 157   0   0
% 1.030000e+02 157   0   0 1.040000e+02 159   0   0
% 1.040000e+02 159   0   0 1.050000e+02 160   0   0
% 1.050000e+02 160   0   0 1.060000e+02 162   0   0
% 1.060000e+02 162   0   0 1.070000e+02 163   0   0
% 1.070000e+02 163   0   0 1.080000e+02 165   0   0
% 1.080000e+02 165   0   0 1.090000e+02 166   0   0
% 1.090000e+02 166   0   0 1.100000e+02 169   0   0
% 1.100000e+02 169   0   0 1.110000e+02 170   0   0
% 1.110000e+02 170   0   0 1.120000e+02 172   0   0
% 1.120000e+02 172   0   0 1.130000e+02 173   0   0
% 1.130000e+02 173   0   0 1.140000e+02 175   1   0
% 1.140000e+02 175   1   0 1.150000e+02 176   3   0
% 1.150000e+02 176   3   0 1.160000e+02 178   5   0
% 1.160000e+02 178   5   0 1.170000e+02 179   7   0
% 1.170000e+02 179   7   0 1.180000e+02 181   9   0
% 1.180000e+02 181   9   0 1.190000e+02 182  11   0
% 1.190000e+02 182  11   0 1.200000e+02 185  15   0
% 1.200000e+02 185  15   0 1.210000e+02 186  17   0
% 1.210000e+02 186  17   0 1.220000e+02 188  18   0
% 1.220000e+02 188  18   0 1.230000e+02 189  20   0
% 1.230000e+02 189  20   0 1.240000e+02 191  22   0
% 1.240000e+02 191  22   0 1.250000e+02 192  24   0
% 1.250000e+02 192  24   0 1.260000e+02 194  26   0
% 1.260000e+02 194  26   0 1.270000e+02 195  28   0
% 1.270000e+02 195  28   0 1.280000e+02 197  30   0
% 1.280000e+02 197  30   0 1.290000e+02 198  32   0
% 1.290000e+02 198  32   0 1.300000e+02 201  35   0
% 1.300000e+02 201  35   0 1.310000e+02 202  37   0
% 1.310000e+02 202  37   0 1.320000e+02 204  39   0
% 1.320000e+02 204  39   0 1.330000e+02 205  41   0
% 1.330000e+02 205  41   0 1.340000e+02 207  43   0
% 1.340000e+02 207  43   0 1.350000e+02 208  45   0
% 1.350000e+02 208  45   0 1.360000e+02 210  47   0
% 1.360000e+02 210  47   0 1.370000e+02 211  49   0
% 1.370000e+02 211  49   0 1.380000e+02 212  51   0
% 1.380000e+02 212  51   0 1.390000e+02 214  52   0
% 1.390000e+02 214  52   0 1.400000e+02 215  54   0
% 1.400000e+02 215  54   0 1.410000e+02 217  56   0
% 1.410000e+02 217  56   0 1.420000e+02 218  58   0
% 1.420000e+02 218  58   0 1.430000e+02 220  60   0
% 1.430000e+02 220  60   0 1.440000e+02 221  62   0
% 1.440000e+02 221  62   0 1.450000e+02 223  64   0
% 1.450000e+02 223  64   0 1.460000e+02 224  66   0
% 1.460000e+02 224  66   0 1.470000e+02 226  68   0
% 1.470000e+02 226  68   0 1.480000e+02 227  69   0
% 1.480000e+02 227  69   0 1.490000e+02 228  71   0
% 1.490000e+02 228  71   0 1.500000e+02 231  75   0
% 1.500000e+02 231  75   0 1.510000e+02 233  77   0
% 1.510000e+02 233  77   0 1.520000e+02 234  79   0
% 1.520000e+02 234  79   0 1.530000e+02 236  81   0
% 1.530000e+02 236  81   0 1.540000e+02 237  83   0
% 1.540000e+02 237  83   0 1.550000e+02 239  85   0
% 1.550000e+02 239  85   0 1.560000e+02 240  86   0
% 1.560000e+02 240  86   0 1.570000e+02 241  88   0
% 1.570000e+02 241  88   0 1.580000e+02 243  90   0
% 1.580000e+02 243  90   0 1.590000e+02 244  92   0
% 1.590000e+02 244  92   0 1.600000e+02 246  94   0
% 1.600000e+02 246  94   0 1.610000e+02 247  96   0
% 1.610000e+02 247  96   0 1.620000e+02 249  98   0
% 1.620000e+02 249  98   0 1.630000e+02 250 100   0
% 1.630000e+02 250 100   0 1.640000e+02 252 102   0
% 1.640000e+02 252 102   0 1.650000e+02 253 103   0
% 1.650000e+02 253 103   0 1.660000e+02 255 105   0
% 1.660000e+02 255 105   0 1.670000e+02 255 107   0
% 1.670000e+02 255 107   0 1.680000e+02 255 109   0
% 1.680000e+02 255 109   0 1.690000e+02 255 111   0
% 1.690000e+02 255 111   0 1.700000e+02 255 115   0
% 1.700000e+02 255 115   0 1.710000e+02 255 117   0
% 1.710000e+02 255 117   0 1.720000e+02 255 119   0
% 1.720000e+02 255 119   0 1.730000e+02 255 120   0
% 1.730000e+02 255 120   0 1.740000e+02 255 122   0
% 1.740000e+02 255 122   0 1.750000e+02 255 124   0
% 1.750000e+02 255 124   0 1.760000e+02 255 126   0
% 1.760000e+02 255 126   0 1.770000e+02 255 128   0
% 1.770000e+02 255 128   0 1.780000e+02 255 130   0
% 1.780000e+02 255 130   0 1.790000e+02 255 132   0
% 1.790000e+02 255 132   0 1.800000e+02 255 136   7
% 1.800000e+02 255 136   7 1.810000e+02 255 137  11
% 1.810000e+02 255 137  11 1.820000e+02 255 139  15
% 1.820000e+02 255 139  15 1.830000e+02 255 141  19
% 1.830000e+02 255 141  19 1.840000e+02 255 143  23
% 1.840000e+02 255 143  23 1.850000e+02 255 145  27
% 1.850000e+02 255 145  27 1.860000e+02 255 147  31
% 1.860000e+02 255 147  31 1.870000e+02 255 149  35
% 1.870000e+02 255 149  35 1.880000e+02 255 151  39
% 1.880000e+02 255 151  39 1.890000e+02 255 153  43
% 1.890000e+02 255 153  43 1.900000e+02 255 156  51
% 1.900000e+02 255 156  51 1.910000e+02 255 158  54
% 1.910000e+02 255 158  54 1.920000e+02 255 160  58
% 1.920000e+02 255 160  58 1.930000e+02 255 162  62
% 1.930000e+02 255 162  62 1.940000e+02 255 164  66
% 1.940000e+02 255 164  66 1.950000e+02 255 166  70
% 1.950000e+02 255 166  70 1.960000e+02 255 168  74
% 1.960000e+02 255 168  74 1.970000e+02 255 170  78
% 1.970000e+02 255 170  78 1.980000e+02 255 171  82
% 1.980000e+02 255 171  82 1.990000e+02 255 173  86
% 1.990000e+02 255 173  86 2.000000e+02 255 175  90
% 2.000000e+02 255 175  90 2.010000e+02 255 177  94
% 2.010000e+02 255 177  94 2.020000e+02 255 179  98
% 2.020000e+02 255 179  98 2.030000e+02 255 181 102
% 2.030000e+02 255 181 102 2.040000e+02 255 183 105
% 2.040000e+02 255 183 105 2.050000e+02 255 185 109
% 2.050000e+02 255 185 109 2.060000e+02 255 187 113
% 2.060000e+02 255 187 113 2.070000e+02 255 188 117
% 2.070000e+02 255 188 117 2.080000e+02 255 190 121
% 2.080000e+02 255 190 121 2.090000e+02 255 192 125
% 2.090000e+02 255 192 125 2.100000e+02 255 196 133
% 2.100000e+02 255 196 133 2.110000e+02 255 198 137
% 2.110000e+02 255 198 137 2.120000e+02 255 200 141
% 2.120000e+02 255 200 141 2.130000e+02 255 202 145
% 2.130000e+02 255 202 145 2.140000e+02 255 204 149
% 2.140000e+02 255 204 149 2.150000e+02 255 205 153
% 2.150000e+02 255 205 153 2.160000e+02 255 207 156
% 2.160000e+02 255 207 156 2.170000e+02 255 209 160
% 2.170000e+02 255 209 160 2.180000e+02 255 211 164
% 2.180000e+02 255 211 164 2.190000e+02 255 213 168
% 2.190000e+02 255 213 168 2.200000e+02 255 215 172
% 2.200000e+02 255 215 172 2.210000e+02 255 217 176
% 2.210000e+02 255 217 176 2.220000e+02 255 219 180
% 2.220000e+02 255 219 180 2.230000e+02 255 221 184
% 2.230000e+02 255 221 184 2.240000e+02 255 222 188
% 2.240000e+02 255 222 188 2.250000e+02 255 224 192
% 2.250000e+02 255 224 192 2.260000e+02 255 226 196
% 2.260000e+02 255 226 196 2.270000e+02 255 228 200
% 2.270000e+02 255 228 200 2.280000e+02 255 230 204
% 2.280000e+02 255 230 204 2.290000e+02 255 232 207
% 2.290000e+02 255 232 207 2.300000e+02 255 236 215
% 2.300000e+02 255 236 215 2.310000e+02 255 238 219
% 2.310000e+02 255 238 219 2.320000e+02 255 239 223
% 2.320000e+02 255 239 223 2.330000e+02 255 241 227
% 2.330000e+02 255 241 227 2.340000e+02 255 243 231
% 2.340000e+02 255 243 231 2.350000e+02 255 245 235
% 2.350000e+02 255 245 235 2.360000e+02 255 247 239
% 2.360000e+02 255 247 239 2.370000e+02 255 249 243
% 2.370000e+02 255 249 243 2.380000e+02 255 251 247
% 2.380000e+02 255 251 247 2.390000e+02 255 253 251];
   
% T = zeros(2*size(A,1),3);
% T(1:2:end,:) = A(:,2:4);
% T(2:2:end,:) = A(:,6:8);

if nargin < 1
    cm_data = cm;
else
    hsv=rgb2hsv(cm);
    hsv(153:end,1)=hsv(153:end,1)+1; % hardcoded
    cm_data=interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,m));
    cm_data(cm_data(:,1)>1,1)=cm_data(cm_data(:,1)>1,1)-1;
    cm_data=hsv2rgb(cm_data);
  
end
end