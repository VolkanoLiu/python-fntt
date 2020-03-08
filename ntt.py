import numpy as np


mod = np.mod


class NttError(Exception):
    def __init__(self, error):
        self.error = error

    def __str__(self, *args, **kwargs):
        return self.error


def ntt(co_array: np.core.multiarray.ndarray, prime: int, primitive_root: int):
    N = len(co_array)
    try:
        divided = mod(prime - 1, N)
        if divided != 0:
            raise(NttError('(prime - 1) show be exactly divided by N'))
    except NttError as err:
        print('\033[1;31m', 'Error: ', err, '\033[0m')
        exit(1)
    value_array = np.zeros(N, dtype=np.uint)
    for k in range(N):
        temp = 0
        for n in range(N):
            temp += mod(co_array[n], prime) * \
                    pow(pow(primitive_root, int((prime - 1) / N)), n * k, prime)
        value_array[k] = mod(temp, prime)
    return value_array


def intt(value_array: np.core.multiarray.ndarray, prime: int, primitive_root: int):
    N = len(value_array)
    try:
        divided = mod(prime - 1, N)
        if divided != 0:
            raise (NttError('(prime - 1) show be exactly divided by N'))
    except NttError as err:
        print('\033[1;31m', 'Error: ', err, '\033[0m')
        exit(1)
    co_array = np.zeros(N, dtype=np.uint)
    for n in range(N):
        temp = 0
        for k in range(N):
            temp += pow(N, prime - 2, prime) * \
                    mod(value_array[k], prime) * \
                    pow(pow(primitive_root, int((prime - 1) / N)), n * k * (prime - 2), prime)
        co_array[n] = mod(temp, prime)
    return co_array


def ntt_core(N: int, prime: int, primitive_root: int):
    try:
        divided = mod(prime - 1, N)
        if divided != 0:
            raise(NttError('(prime - 1) show be exactly divided by N'))
    except NttError as err:
        print('\033[1;31m', 'Error: ', err, '\033[0m')
        exit(1)
    core = np.zeros((N, N), dtype=np.int)
    for k in range(N):
        temp = 0
        for n in range(N):
            core[k][n] = pow(pow(primitive_root, int((prime - 1) / N)), n * k, prime)
    return core


def ntt_tr_core(N: int, prime: int, primitive_root: int):
    try:
        divided = mod(prime - 1, N)
        if divided != 0:
            raise(NttError('(prime - 1) show be exactly divided by N'))
    except NttError as err:
        print('\033[1;31m', 'Error: ', err, '\033[0m')
        exit(1)
    core = np.zeros((N, N), dtype=np.int)
    phi = pow(primitive_root, int((prime - 1) / N / 2), prime)
    for k in range(N):
        temp = 0
        for n in range(N):
            core[k][n] = pow(pow(primitive_root, int((prime - 1) / N)), n * k, prime) * pow(phi, n)
    return mod(core, prime)


def intt_core(N: int, prime: int, primitive_root: int):
    try:
        divided = mod(prime - 1, N)
        if divided != 0:
            raise (NttError('(prime - 1) show be exactly divided by N'))
    except NttError as err:
        print('\033[1;31m', 'Error: ', err, '\033[0m')
        exit(1)
    core = np.zeros((N, N), dtype=np.int)
    for n in range(N):
        temp = 0
        for k in range(N):
            core[k][n] = pow(N, prime - 2, prime) * \
                    pow(pow(primitive_root, int((prime - 1) / N)), n * k * (prime - 2), prime)
    return np.mod(core, prime)


# ntt_arr = np.arange(1, 9, dtype=np.uint)
# print(np.polydiv(np.convolve([1, 2, 3, 4], [1, 2, 3, 4]), [1, 0, 0, 0, 1]))
# # print(ntt_arr)
# print(ntt(ntt_arr, 17, 3))
# NN = len(ntt_arr)
# phi = pow(3, int((17 - 1) / NN / 2), 17)
# for n in range(NN):
#     ntt_arr[n] = mod(ntt_arr[n] * pow(phi, n), 17)
# # print(ntt_arr)
# print(ntt(ntt_arr, 17, 3))
#
#
# exit(0)
