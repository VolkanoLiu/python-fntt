import fntt
import ntt
import numpy as np
import matplotlib.pyplot as plt

# test_array_ntt = np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=int)
# test_array_fntt = np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int64)
# for i in range(32):
#     test_array_fntt = np.random.randint(0, 65537, int(65536 / 2), dtype=np.int64)
#     # print(test_array_fntt)
#     b_test = fntt.fntt(test_array_fntt, 65537, 3)
#     # print(b_test)
#     res = fntt.fntt(b_test, 65537, 3, True)
#     print(np.all(test_array_fntt == res))


test_array = np.arange(1, 9, dtype=np.int64)
# test_array = np.random.randint(0, 65537, 1024, dtype=np.int64)
a = fntt.linear_conv(test_array, test_array, 17, 3)

