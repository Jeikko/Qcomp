#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""
Tools used to manipulate continued fractions

Function names are converters from one number representation to another.
bits: list of bits, the most significant one coming first.
dec: float.
cfrac: list of integers appearing in a continued fraction. The first element is
the first integer, [a, b, c] corresponds to 1(a+1/(b+1/(c))).
frac: a pair (numerator, denominator)
"""

def bits_to_dec(l):
    """
    Convert a number from a bit to decimal (float) representation.
    
    l is a list of bits, the most significant one coming first.
    bits_to_dec([]) -> 0
    """
    if l == []:
        return 0
    else:
        return 0.5*l[0] + 0.5*bits_to_dec(l[1:])

def dec_to_cfrac(d, t):
    """
    Give the continued fraction decomposition of a float.
    
    d is the number to decompose, 0<t<1 is a threshold. If a number lower than
    this value appears in the decomposition, the algorithm truncates it and
    stops. 
    """
    if d < t:
        return []
    return [int(1/d)] + dec_to_cfrac(1/d-int(1/d), t)

def cfrac_to_frac(l):
    """
    Give the numerator and denominatior of a continued fraction.
    
    The function returns a pair (numerator, denominatior)
    """
    num = 0
    den = 1
    for r in reversed(l):
        num, den = den, num + r*den
    return num, den

def dec_to_frac(d, t):
    """
    Turn a decimal number into a fraction
    
    The continued fraction decomposition (with a truncation) is used.
    d is the number to decompose, 0<t<1 is a threshold. If a number lower than
    this value appears in the decomposition, the algorithm truncates it and
    stops.
    """
    return cfrac_to_frac(dec_to_cfrac(d, t))

def bits_to_frac(l, t):
    """
    Turn a binary number into a fraction
    
    The continued fraction decomposition (with a truncation) is used.
    b is the number to decompose, a list of bits (the most significant one
    coming first). 0<t<1 is a threshold. If a number lower than this value
    appears in the decomposition, the algorithm truncates it and stops.
    """
    return dec_to_frac(bits_to_dec(l), t)
