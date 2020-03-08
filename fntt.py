#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    fntt.py
# @Author:      volkanoLiu
# @Time:        2020-03-08
# @License:

import numpy as np
import matplotlib.pyplot as plt


def print_err_exit(err: str):
    print('\033[1;31m', 'Error: ', err, '\033[0m')
    exit(1)


def get_2_power(num: int):
    if num & (num - 1):
        print_err_exit('length is not 2^n')

    digit = 0
    num_temp = num
    while num_temp != 1:
        num_temp = num_temp >> 1
        digit += 1
    return digit


def fntt_core(point_array: np.ndarray, prime: int, primitive_root: int, index_array: np.ndarray, inverse=False):
    length = len(index_array)
    if length == 1:
        return point_array[index_array[0]]
    else:
        multiply_serial = np.zeros(int(length / 2), dtype=np.int64)
        if not inverse:
            for i in range(int(length / 2)):
                multiply_serial[i] = pow(primitive_root, int(i * (prime - 1) / length), prime)
        else:
            for i in range(int(length / 2)):
                # primitive_root_i = pow(primitive_root, prime - 2, prime)
                multiply_serial[i] = pow(primitive_root, int(i * (prime - 1) * (prime - 2) / length), prime)

        even_half = index_array[0::2]
        odd_half = index_array[1::2]
        even_half_fntt = fntt_core(point_array, prime, primitive_root, even_half, inverse)
        odd_half_fntt = fntt_core(point_array, prime, primitive_root, odd_half, inverse)
        head_half = np.mod(even_half_fntt + multiply_serial * odd_half_fntt, prime)
        tail_half = np.mod(even_half_fntt - multiply_serial * odd_half_fntt, prime)

    return np.append(head_half, tail_half)


def fntt(input_array: np.ndarray, prime: int, primitive_root: int, inverse=False) -> np.ndarray:
    """Fast Number Theory Transform

    :param input_array: input array
    :param prime: a prime number
    :param primitive_root: the primitive root of the prime
    :param inverse: if True, take the IFNTT
    :return: the FNTT/IFNTT result
    """
    if input_array.dtype != np.int64:
        print_err_exit('invalid data type')
    data_length = len(input_array)
    get_2_power(data_length)

    index_array = np.arange(0, data_length, dtype=np.int64)
    result = fntt_core(input_array, prime, primitive_root, index_array, inverse)
    if inverse:
        result = np.mod(result * pow(data_length, prime - 2, prime), prime)

    return result


def linear_conv(input_array_a: np.ndarray, input_array_b: np.ndarray, prime: int, primitive_root: int):
    """calculate the convolution of two arrays

    :param input_array_a:
    :param input_array_b:
    :param prime:
    :param primitive_root:
    :return:
    """
    len_a = len(input_array_a)
    len_b = len(input_array_b)
    if len_a != len_b:
        print_err_exit("not the same length!")
    get_2_power(len_a)
    phi = pow(primitive_root, int((prime - 1) / len_a / 2), prime)
    a_tr = np.zeros(len_a, dtype=np.int64)
    b_tr = np.zeros(len_b, dtype=np.int64)
    for n in range(len_a):
        a_tr[n] = np.mod(input_array_a[n] * pow(phi, n, prime), prime)
        b_tr[n] = np.mod(input_array_b[n] * pow(phi, n, prime), prime)
    a_tr_fntt = fntt(a_tr, prime, primitive_root)
    b_tr_fntt = fntt(b_tr, prime, primitive_root)
    a_multiply_b_tr_fntt = a_tr_fntt * b_tr_fntt
    a_multiply_b_tr = fntt(a_multiply_b_tr_fntt, prime, primitive_root, True)
    phi_reverse = pow(primitive_root, int((prime - 1) / len_a / 2 * (prime - 2)), prime)
    a_multiply_b = np.zeros(len_a, dtype=np.int64)
    for n in range(len_a):
        a_multiply_b[n] = np.mod(a_multiply_b_tr[n] * pow(phi_reverse, n, prime), prime)
    return a_multiply_b
