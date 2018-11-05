#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""
Quantum computer simulation v1

The purpose of this file is to check my correct understanding of the bases of 
quantum computer simulation.

The simulation is based on a single class, Qreg (quantum register), which
represents a quantum register made of N qubits.
Several methods acting on that register (gates, measurements, but also direct
and reverse QFT and adders) are defined.

That class is then used to implement the phase estimation technique. The
quantum part of Shor algorithm is also implemented and applied.
"""


import numpy
import math
import random

#import binpad, get_bit, add_zeroes and bit_spread
from bit_ops import *
from cont_frac import bits_to_dec, dec_to_frac 

def cExp(x):
    """Return exp(i x)"""
    return math.cos(x) + 1.J * math.sin(x)

#===============================================================================
# Standard gate matrices
# Ud is U dagger, the adjoint of U
#===============================================================================
H = numpy.array([[1,1],[1,-1]])/math.sqrt(2)
X = numpy.array([[0,1],[1,0]])
Y = numpy.array([[0,-1.J],[1.J,0]])
Z = numpy.array([[1,0],[0,-1]])
S = numpy.array([[1,0],[0,1.J]])
Sd = numpy.array([[1,0],[0,-1.J]])
T = numpy.array([[1,0],[0,(1.+1.J)/math.sqrt(2)]])
Td = numpy.array([[1,0],[0,(1.-1.J)/math.sqrt(2)]])
# These matrices are the square root of X, used for Toffoli gates
srX = numpy.array([[(1-1.j)/2, (1+1.j)/2], [(1+1.j)/2, (1-1.j)/2]])
srXd = numpy.array([[(1+1.j)/2, (1-1.j)/2], [(1-1.j)/2, (1+1.j)/2]])


class Qreg:
    """
    Quantum register class
    
    r = Qreg(N, init=b) creates a quantum register of N qubits (that is, a
    list of 2^N complex amplitudes) and puts all qubits in the state b (default
    value 0).
    
    r.gate and r.cgate (controlled gate) are general methods that apply a gate
    to some qubits of the register. r.swap swaps two qubits.
    r.X, r.H, r.cnot and r.Toffoli are shorthands for these commands.
    r.QFT and r.QFTd apply direct and reverse QFT to some bits of the register
    
    r.measurement performs a measurement on a given qubit, modifies the state
    of the register according to the result and returns that result. To keep
    successive measurements consistent for the user, a list of all measured
    qubits is kept. Without it, measuring all qubits of a 2-qubit register
    would be r.measurement(0) and r.measurement(0) again (after the first
    measurement, the register only has one qubit)
    
    It is also possible to append a second register: r.append(s) is a register
    of N + M qubits if s has M qubits.
    
    In the print() representation of the states, the highest-weight qubit is
    on the left (as in the binary representation of an integer)
    """
    def __init__(self, N, init=0):
        """
        Create a quantum register of N qubits
        
        Initialize all of them in the state init
        """
        self.N = N
        self.amplitudes = [0. for _ in range(2**N)]
        if init == 0:
            self.amplitudes[0] = 1.
        if init == 1:
            self.amplitudes[-1] = 1.
        #self.aq: active qubits, i.e. those that have not been measured
        self.aq = list(range(self.N))
        #self.dq: dropped (= measured) qubits
        self.dq = []
    def _cdq(self, i):
        """
        Convert qubit indices, taking into account dropped qubits
        
        If the qubit 0 has been dropped, cdq will convert the qubit i (as
        seen by the user) to i-1 (as stored in the amplitude list).
        That method is already taken into account in the gate/swap methods,
        so the user probably doesn't need it.
        i can be a qubit index, and the method will return the modified index
        It can also be a list of indices, and the method will return the list
        of modified indices
        
        Raise a ValueError if i is a qubit that has already been dropped
        """
        if type(i) == list:
            return [self._cdq(j) for j in i]
        count = 0
        for j in self.dq:
            if i == j:
                raise ValueError("The qubit {} has been dropped".format(i))
            if j < i:
                count += 1
        return i - count
    def __str__(self):
        """
        Display the relevant states of a register
        """
        l_amp = ["{} |{}> (sqm {})".format(self.amplitudes[i],
                                           binpad(i, self.N),
                                           abs(self.amplitudes[i])**2)
                  for i in range(2**self.N)
                  if abs(self.amplitudes[i])**2 > 10**(-4)]
        if len(self.dq) != 0:
            pdq = "\nDropped qubits: " + ", ".join([str(i) for i in self.dq])
        else:
            pdq = ""
        return "State of {} active and {} dropped qubits:\n".format(
                self.N, len(self.dq)
                ) + "\n".join(l_amp) + pdq + '\n'
    def append(self, other):
        """
        Append another register
        
        other will not be modified, but its dropped qubits won't be transferred
        to the current register
        """
        temp = [elem for elem in self.amplitudes]
        #Extend the list of amplitudes
        self.amplitudes = [0. for _ in range(2**(self.N + other.N))]
        for j in range(2**other.N):
            for i in range(2**self.N):
                self.amplitudes[j * 2**self.N + i] = \
                other.amplitudes[j] * temp[i]
        #Add active qubits. The dropped qubits are not modified
        self.aq += list(range(
            self.N + len(self.dq), 
            self.N + len(self.dq) + other.N))
        self.N += other.N 
    
    #===========================================================================
    # General gate methods
    #===========================================================================
    def gate(self, op, T):
        """
        Apply a gate to target qubits
        
        op is a quare matrix of size 2^Lx2^L
        T is a list of L indices, or can be a single index if op is 2x2
        """
        T = self._cdq(T)
        if type(T) != list:
            T = [T]
#         Instead of creating a 2^Nx2^N matrix for op, divide the problem into
#         the computation of 2^(N-L) matrix products, for each of the 2^(N-L)
#         possible values for the unmodified (not in T) qubits
        for i in range(2**(self.N-len(T))):
#             mask is the number whose binary representation is i on all bits
#             except the ones in T, and 0 for them
            mask = add_zeroes(i, T)
#             bit_spread(j, T) is the number whose binary representation is 0
#             on all bits except the ones in T, and the bits of j for them.
#             J is the list of indices where the operator op should apply
            J = [mask + bit_spread(j, T)
                 for j in range(2**len(T))]
            amp = numpy.array([[self.amplitudes[elem]] for elem in J])
#             Modified amplitudes
            res = op.dot(amp)
#             Updating the register
            for j in range(2**len(T)):
                self.amplitudes[J[j]] = res[j][0]
    def cgate(self, op, control, T):
        """
        Apply a gate to target qubits, controlled by another one
        
        op is a quare matrix of size 2^Lx2^L
        control is the index of the controlling qubit
        T is a list of L indices, or can be a single index if op is 2x2
        """
        control = self._cdq(control)
        T = self._cdq(T)
        if type(T) != list:
            T = [T]
#         Instead of creating a 2^Nx2^N matrix for op, divide the problem into
#         the computation of 2^(N-L-1) matrix products, for each of the
#         2^(N-L-1) possible values for the unmodified (not in T, not control)
#         qubits
        for i in range(2**(self.N - len(T) - 1)):
#             mask is the number whose binary representation is i on all bits
#             except the ones in T and control, and 0 for them
            mask = add_zeroes(i, T + [control])
#             bit_spread(j, T) is the number whose binary representation is 0
#             on all bits except the ones in T, and the bits of j for them.
#             J is the list of indices where the operator op should apply.
#             The control bit is put to 1, because the gate will not modify the
#             states where it is 0
            J = [mask + 2**control + bit_spread(j, T)
                 for j in range(2**len(T))]
            amp = numpy.array([[self.amplitudes[elem]] for elem in J])
#             Modified amplitudes
            res = op.dot(amp)
#             Updating the register
            for j in range(2**len(T)):
                self.amplitudes[J[j]] = res[j][0]
    def swap(self, t1, t2):
        """Swap the qubits t1 and t2"""
        t1 = self._cdq(t1)
        t2 = self._cdq(t2)
#         Very similar to r.gate, because swapping is applying the swap
#         operation to them.It is simply written explicitely. 
        for i in range(2**(self.N-2)):
            mask = add_zeroes(i, [t1, t2])
            j1 = mask + bit_spread(1, [t1, t2])
            j2 = mask + bit_spread(2, [t1, t2])
            temp = self.amplitudes[j1]
            self.amplitudes[j1] = self.amplitudes[j2]
            self.amplitudes[j2] = temp
    
    #===========================================================================
    # Specific gates
    #===========================================================================
    def cnot(self, control, target):
        """Apply the controlled-not gate"""
        self.cgate(X, control, [target])
    def H(self, target):
        """Apply the Hadamard gate"""
        self.gate(H, [target])
    def X(self, target):
        """Apply the X/not gate"""
        self.gate(X, [target])
    def Toffoli(self, c1, c2, t):
        """Apply the Toffoli gate"""
        self.cgate(srX, c2, t)
        self.cnot(c1, c2)
        self.cgate(srXd, c2, t)
        self.cnot(c1, c2)
        self.cgate(srX, c1, t)
    def half_adder(self, a, b, s, c):
        """
        Half-adder operation
        
        s and c are two result qubits, initialized to 0.
        s is put to a XOR b, c to a AND b
        """
        self.Toffoli(a, b, c)
        self.cnot(a, s)
        self.cnot(b, s)
    def QFT(self, l_qb=None):
        """
        Apply the quantum Fourier transform to some qubits
        
        l_qb is the list of qubits where the transform is to be applied.
        Default value is all qubits
        It is done by repeatedly applying to all qubits first a Hadamard gate,
        then several rotations around z. Qubits are swapped at the end
        """
        def Rk(k):
            return numpy.array([[1, 0], [0, cExp(2 * math.pi / 2**k)]])
        if l_qb == None:
            l_qb = self.aq
        for i in range(len(l_qb)):
            t = l_qb[i]
            self.H(t)
            for j in range(i+1, len(l_qb)):
                c = l_qb[j]
                self.cgate(Rk(j-i+1), c, t)
        for i in range(len(l_qb)//2):
            self.swap(l_qb[i], l_qb[-i-1])
    def QFTd(self, l_qb=None):
        """
        Apply the reverse quantum Fourier transform to some qubits
        """
        def Rkd(k):
            return numpy.array([[1, 0], [0, cExp(-2 * math.pi / 2**k)]])
        if l_qb == None:
            l_qb = self.aq
        for i in range(len(l_qb)//2):
            self.swap(l_qb[i], l_qb[-i-1])
        for i in reversed(range(len(l_qb))):
            t = l_qb[i]
            for j in reversed(range(i+1, len(l_qb))):
                c = l_qb[j]
                self.cgate(Rkd(j-i+1), c, t)
            self.H(t)
        
    
    #===========================================================================
    # Measurement methods
    #===========================================================================
    def measurement(self, target):
        """
        Perform a projective measurement on a given qubit
        
        The measurement is performed in the basis |0>, |1>. target is the
        qubit to be measured. The result is random and is the return value of
        the method.
        The method also updates the quantum register in accordance to the
        measurement result.
        """
        target_c = self._cdq(target)
#         w0 and w1 are the cumulative weights associated to the two possible
#         results. At the end of the for loop, w0 + w1 = 1
        w0 = 0
        w1 = 0
        for i in range(2**self.N):
            if get_bit(i, target_c) == 0:
                w0 += abs(self.amplitudes[i])**2
            else:
                w1 += abs(self.amplitudes[i])**2
#         Pick a random result, according to the weights
        if random.random() < w0:
            res = 0
            w = w0
        else:
            res = 1
            w = w1
#         Update the amplitudes
        for i in range(2**(self.N - 1)):
            self.amplitudes[i] = self.amplitudes[
                add_zeroes(i, target_c) + res * 2**target_c]/math.sqrt(w)
#         Update other Qreg values
        self.N -= 1
        self.dq.append(target)
        self.aq.remove(target)
#         Clean the amplitude list
        del self.amplitudes[2**self.N:]
        return res
    def s_measurements(self, T=None):
        """
        Perform successive measurements and return the result
        
        If T is not given, measurements are done on all qubits.
        """
        if T == None:
            T = [elem for elem in self.aq]
        return [self.measurement(t) for t in T]


def phase_estimation(U, r_u, n):
    """
    Implementation of the phase estimation circuit
    
    U is the (2^Nx2^N) operator whose eigenvalue is to be measured
    r_u is an auxiliary register prepared in the correct eigenstate for u. It
    has N qubits
    n is the number of qubits that will be used in the first register
    
    Return the first register, after a measurement has been performed on the
    auxiliary one to discard it.
    
    Note that r_u is copied but is not modified by the function.
    """
    aux_qb = list(range(n, n + r_u.N))
    def Upow(j):
        return numpy.linalg.matrix_power(U, 2**j)
    r = Qreg(n)
    for i in range(n):
        r.H(i)
    r.append(r_u)
    for i in reversed(range(n)):
        r.cgate(Upow(n-i-1), i, aux_qb)
#     Principle of implicit measurement allows to measure the auxiliary qbits
#     here to reduce the system size and give a faster QFT
    r.s_measurements(aux_qb)
    r.QFTd(r.aq)
    return r

#===============================================================================
# Example of phase estimation
#===============================================================================
def example_phase_estimation():
    """Example of phase estimation"""
    U = numpy.array([[1, 0], [0, cExp(2 * math.pi * 0.3)]])
    d = Qreg(1, init=1)
    r = phase_estimation(U, d, 10)
    res = r.s_measurements()
    print(bits_to_dec(res))

#===============================================================================
# Implementation of Shor's algorithm
#===============================================================================
def Shor(N=21, a=10):
    """
    Implement the quantum part of Shor algorithm.
    
    Note that U is built as a matrix, line-by-line, instead of being a series
    of gate operations. It's simpler and still allows to understand and
    implement the algorithm, but it loses the O(L^3) efficiency. Implementing
    an actual construction for U requires better tools to work on circuits,
    which is beyond the scope of this version
    """
    L = int(math.ceil(math.log(N, 2)))
    
    # Instead of building it with arithmetic operations, U is seen and built as
    # a matrix here. It's simpler and still allows to understand and implement
    # the algorithm, but it loses the O(L^3) efficiency
    U = numpy.zeros([2**L, 2**L])
    for i in range(N):
        U[a*i%N][i] = 1
    for i in range(N, 2**L):
        U[i][i] = 1
    
    # Preparing the auxiliary register
    d = Qreg(L)
    d.X(0)
    
    # Running the algorithm a number of times, giving each time the fraction
    # approximation of s/r. The order should be the least common multiple of
    # all denominators 
    for i in range(10):
        r = phase_estimation(U, d, 10)
        res = bits_to_dec(r.s_measurements())
        print(dec_to_frac(res, 0.1), res)

Shor(15, 7)

