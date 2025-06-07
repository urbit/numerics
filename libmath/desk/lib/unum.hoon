/+  *math,  :: notably, this shadows floating-point doors
    twoc
  ::
::::  Mathematical library
::
::  UNIVERSAL NUMBERS (UNUMS)
::
::  Type III universal numbers (unums) are divided into three classes:
::  1. Posits, or precise real numbers (well, as real as floats anyway).
::  2. Quires, or 
::  3. Valids, or bounded intervals.
::
::  We represent these using Hoon auras at three bitwidths:
::  - @rpb, @rph, @rps  :: posits
::  - @rqb, @rqh, @rqs  :: quires
::  - @rvb, @rvh, @rvs  :: valids
::
::  The bitwidths are 8-bit Bytes, 16-bit Halfs, and 32-bit Singles.
::  While we square on quires and valids, the current implementation is
::  only for posits.
::
::  - Gustafson & Yonemoto (2017), "Beating Floating Point at its Own Game:
::    Posit Arithmetic", Supercomputing Frontiers and Innovations. 4 (2).
::    Publishing Center of South Ural State University, Chelyabinsk, Russia.
::    http://www.johngustafson.net/pdfs/BeatingFloatingPoint.pdf
|%
::  Representation of a Type III Unum Posit
::  (uq, uv)
+$  up
  $%  $:  %p    :: real-valued posit
          b=@u  :: bitwidth (bloq), ∈ 3 (byte), 4 (half), 5 (single)
          s=?   :: sign, 0 (+) or 1 (-)
          r=@s  :: regime
          e=@s  :: exponent
          f=@u  :: fraction
      ==
      [%n ~]    :: Not a Real (NaR), unum NaN
  ==
::
+$  pb
  $%  s=?         :: sign bit
      r=[w=@ p=@] :: regime bits, bitwidth first
      e=[w=@ p=@] :: exponent bits, bitwidth first
      f=@         :: fraction bits, remaining to total bitwidth
  ==
::  Type III Unum Posit, 8-bit width ("byte")
++  rpb
  ^|
  |_  $:  rtol=_.1e-5         :: relative tolerance for precision of operations
      ==
  |%
  ::
  ::  Operations
  ::
  ::    +sea:  @pb -> $up
  ::
  ::  Returns the +$up representation of a floating-point atom.
  ::    Examples
  ::      :: posit 1.0
  ::      > (sea 0b100.0000)
  ::      [%p s=%.y r=0 e=0 f=1]
  ::      :: posit 0.5
  ::      > (sea 0b11.1000)
  ::      [%p s=%.y r= e= f=0b0]
  ::      :: posit largest possible 8-bit value, 2**24
  ::      > (sea 0b111.1111)
  ::      [%p s=%.y r=0b111.1111 e=0b0 f=0b0]
  ::      :: posit largest possible negative 8-bit value, -2**24+1
  ::  Source
  ++  sea  !!
  ::
  ::    +to-rs:  @pb -> @rs
  ::
  ::  Returns the closest @rs equivalent to a posit.
  ::  $((1-3s)+f)\times 2^{(1-2s)\times(4r+e+s)}$
  ::
  ::  Source
  ++  to-rs
    |=  =@rpb
    ^-  @rs
    =/  =pb  (sea rpb)
    =/  s=@  `@`s.pb
    =/  r=@  (get-regime r.pb)
    =/  e=@  (sun:^rs (to-sd:twoc e.pb))
    =/  f=@  (sun:^rs (to-sd:twoc f.pb))
    %+  mul:rs
      %+  add:rs
        (sub:rs .1 (mul:rs .3 s))
      f
    %-  pow-n:rs
    %+  mul:rs
      :(sub:rs .1 s s)
    :(add:rs (mul:rs .4 r) e s)
  --
::  Type III Unum Posit, 16-bit width ("half")
::  Type III Unum Posit, 32-bit width ("single")
--
