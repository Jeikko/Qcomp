#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""
Offer some functions to manipulate numbers, bit-wise

binpad -> display a number in binary form
get_bit -> return the n-th bit of a number
add_zeroes -> insert zeroes in the binary representation of a number
bit_spread -> create a new number by moving the bits of another number
"""

def binpad(i, L):
    """
    Return the bit representation of a number in a given length.
    
    i is the number to be written, L is the length of the result.
    The result is a string of '0' and '1'.
    Raise an AssertionError if L is too short for i.
    """
    assert i < 2**L
    return bin(i).split("b")[1].rjust(L, "0")

def get_bit(i, n):
    """
    Return the nth bit of number i
    """
    return (i >> n)%2

def add_zeroes(i, lpos):
    """
    Insert zeroes in the binary representation of a number.
    
    i is the number to be altered, lpos is the position(s) where zeroes
    should be added, either an int or a list of int.
    Positions are counted from the right, starting from 0.
    The result is an int.
    
    Example:
    add_zeroes(15, 1) will turn 15 (0b1111) into 29 (0b11101)
    add_zeroes(15, [1, 3]) will turn 15 into 53 (0b110101)
    """
    if type(lpos) != list:
        lpos = [lpos]
    for pos in sorted(lpos):
        q, r = i//2**pos, i%2**pos
        i = q*2**(pos+1) + r
    return i

def bit_spread(j, lpos):
    """
    Create a new number by moving the bits of another number.
    
    i is the number whose bits are to be spread, lpos is the position,
    or a list of positions, where the bits should be put.
    The bits whose position is not in lpos are set to 0.
    The bit of least weight is sent to the first element in lpos, and so
    on.
    Positions are counted from the right, starting from 0.
    The result is an int.
    Raise an AssertionError if there are not enough bit positions
    
    Example:
    bit_spread(0, [3, 1]) (0b11) ->  0 (0b0000)
    bit_spread(1, [3, 1]) (0b11) ->  8 (0b1000)
    bit_spread(2, [3, 1]) (0b11) ->  2 (0b0010)
    bit_spread(3, [3, 1]) (0b11) -> 10 (0b1010)
    bit_spread(4, [3, 1]) (0b11) -> AssertionError (4 is 3 bits long)
    """
    if type(lpos) != list:
        lpos = [lpos]
    assert j < 2**len(lpos)
    nb = 0
    for pos in lpos:
        nb += j%2 * 2**pos
        j = j//2
    return nb
