#!/usr/bin/env python3
"""Oracle for /lib/complex and Lagoon %cplx, against NumPy complex64/complex128.

A complex value is packed as a single atom with the REAL component in the low
2^(bloq-1) bits and the IMAGINARY component in the high bits -- matching
SoftBLAS complexN_t{real;imag} over little-endian bytes.  @cs (bloq 6) = two
float32; @cd (bloq 7) = two float64.

Prints the packed hex the Hoon tests assert.  Run:  python3 complex_check.py
"""
import numpy as np, struct


def pack(z, n):
    """complex -> packed @c hex; n = component width in bits (32 or 64)."""
    z = (np.complex64 if n == 32 else np.complex128)(z)
    if n == 32:
        re = struct.unpack('<I', struct.pack('<f', float(z.real)))[0]
        im = struct.unpack('<I', struct.pack('<f', float(z.imag)))[0]
    else:
        re = struct.unpack('<Q', struct.pack('<d', float(z.real)))[0]
        im = struct.unpack('<Q', struct.pack('<d', float(z.imag)))[0]
    return (im << n) | re


def grp(v):
    """atom -> dotted hex (4-hex-digit groups), Hoon style."""
    s = '%x' % v
    out = ''
    while len(s) > 4:
        out = '.' + s[-4:] + out
        s = s[:-4]
    return '0x' + s + out


def show(label, z, n=32):
    print('  %-26s %s' % (label, grp(pack(z, n))))


def main():
    print('@cs (complex-single, float32 components):')
    a, b = np.complex64(1 + 2j), np.complex64(3 + 4j)
    show('1+2i', 1 + 2j)
    show('3+4i', 3 + 4j)
    show('add (1+2i)+(3+4i)=4+6i', a + b)
    show('sub =-2-2i', a - b)
    show('mul =-5+10i', a * b)
    show('div (2+2i)/(1+1i)=2', np.complex64(2 + 2j) / np.complex64(1 + 1j))
    show('conj(1+2i)=1-2i', np.conj(a))
    show('abs(3+4i)=5', np.complex64(abs(b)))
    show('one 1+0i', 1 + 0j)
    show('cumsum [1+2i,3+4i,-2-2i]=2+4i', a + b + np.complex64(-2 - 2j))
    x = np.array([1 + 2j, 3 + 4j], dtype=np.complex64)
    y = np.array([5 + 6j, 7 + 8j], dtype=np.complex64)
    show('dot (bilinear) sum x*y', np.complex64(np.sum(x * y)))
    show('dotc (hermitian) vdot(x,y)', np.complex64(np.vdot(x, y)))
    print('  (mmul: [[1+1i,2][0,1]]x[[1,0][0,1+1i]] = [[1+1i,2+2i][0,1+1i]])')

    print('@cd (complex-double, float64 components):')
    show('1+2i', 1 + 2j, 64)
    show('3+4i', 3 + 4j, 64)
    show('add 4+6i', 4 + 6j, 64)
    show('mul -5+10i', np.complex128(1 + 2j) * np.complex128(3 + 4j), 64)
    show('div=2', np.complex128(2 + 2j) / np.complex128(1 + 1j), 64)
    show('one 1+0i', 1 + 0j, 64)


if __name__ == '__main__':
    main()
