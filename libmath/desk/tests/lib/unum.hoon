  ::  /tests/lib/unum
::::
::    Posits
::
/+  math,
    *test,
    unum
^|
|%
::  posit8 tests
++  test-values-rpb  ^-  tang
  ;:  weld
    ::  0
    %+  expect-eq
      !>  `@rpb`0b0
      !>  zero:rpb:unum
    ::  1
    %+  expect-eq
      !>  `@rpb`0b100.0000
      !>  one:rpb:unum
    ::  -1
    %+  expect-eq
      !>  `@rpb`0b1100.0000
      !>  neg-one:rpb:unum
    ::  NaR
    %+  expect-eq
      !>  `@rpb`0b1000.0000
      !>  nar:rpb:unum
    ::  pi
    %+  expect-eq
      !>  `@rpb`0b110.1001
      !>  pi:rpb:unum
    ::  tau
    %+  expect-eq
      !>  `@rpb`0b111.0101
      !>  tau:rpb:unum
    ::  e
    %+  expect-eq
      !>  `@rpb`0b110.0110
      !>  e:rpb:unum
    ::  phi
    %+  expect-eq
      !>  `@rpb`0b101.0100
      !>  phi:rpb:unum
    ::  sqrt(2)
    %+  expect-eq
      !>  `@rpb`0b100.1101
      !>  sqt2:rpb:unum
    ::  TODO other constants
    ::  huge
    %+  expect-eq
      !>  `@rpb`0b111.1111
      !>  huge:rpb:unum
    ::  neg-huge
    %+  expect-eq
      !>  `@rpb`0b1000.0001
      !>  neg-huge:rpb:unum
    ::  tiny
    %+  expect-eq
      !>  `@rpb`0b1
      !>  tiny:rpb:unum
    ::  neg-tiny
    %+  expect-eq
      !>  `@rpb`0b1111.1111
      !>  neg-tiny:rpb:unum
  ==
::
++  test-from-rpb  ^-  tang
  ;:  weld
    ::  0
    %+  expect-eq
      !>  [%z 3 ~]
      !>  (from:rpb:unum zero:rpb:unum)
    ::  1
    %+  expect-eq
      !>  [%p 3 s=%.y r=--0 e=--0 f=0]
      !>  (from:rpb:unum one:rpb:unum)
  ==
::
++  test-into-rpb  ^-  tang
  ;:  weld
    ::  0
    %+  expect-eq
      !>  zero:rpb:unum
      !>  (into:rpb:unum [%z 3 ~])
    ::  1
    %+  expect-eq
      !>  one:rpb:unum
      !>  (into:rpb:unum [%p 3 s=%.y r=--0 e=--0 f=0])
    ::  1.25
    %+  expect-eq
      !>  `@rpb`0b100.1000
      !>  (into:rpb:unum [%p 3 s=%.y r=--0 e=--0 f=0b1000])
    ::  0.5
    %+  expect-eq
      !>  `@rpb`0b10.0000
      !>  (into:rpb:unum [%p 3 s=%.y r=-1 e=--0 f=0])
    ::  -1.0
    %+  expect-eq
      !>  neg-one:rpb:unum
      !>  (into:rpb:unum [%p 3 s=%.n r=--0 e=--0 f=0])
    ::  huge
    %+  expect-eq
      !>  huge:rpb:unum
      !>  (into:rpb:unum [%p 3 s=%.y r=--6 e=--0 f=0b0])
    ::  neg-huge
    %+  expect-eq
      !>  neg-huge:rpb:unum
      !>  (into:rpb:unum [%p 3 s=%.n r=-6 e=--0 f=0b0])
    ::  tiny
    %+  expect-eq
      !>  tiny:rpb:unum
      !>  (into:rpb:unum [%p 3 s=%.y r=-6 e=--0 f=0b0])
    ::  neg-tiny
    %+  expect-eq
      !>  neg-tiny:rpb:unum
      !>  (into:rpb:unum [%p 3 s=%.n r=--6 e=--0 f=0b0])
  ==
::
::  round-trip tests

:: 1.25

:: s = 0
:: r = 0  <-  k = 1, R_0 = 1
:: e = 0
:: f = 01000

:: ((1 - 3*0) + 2**(-5)*(1*2**3)) * 2**((1-2*0)*(1*0+0+0))
::     1      + (1/32) * 8        * 1
::     1      + 1/4               * 1

:: ---

:: 0.5

:: s = 0
:: r = -1 <-  k = 1, R_0 = 0
:: e = 0
:: f = 0

:: ((1 - 3*0) + 2**(-5)*(0*2**3)) * 2**((1-2*0)*(1*-1+0+0))

:: ---

:: 1.0

:: s = 0
:: r = 0  <- k = 1, R_0 = 1
:: e = 0
:: f = 0

:: ((1 - 3*0) + 2**(-5)*(0*2**3)) * 2**((1-2*0)*(1*0+0+0))

:: ---

:: -1.0

:: s = 1
:: r = 0  <- k = 1, R_0 = 1
:: e = 0
:: f = 0

:: ((1 - 3*1) + 2**(-5)*(0*2**3)) * 2**((1-2*1)*(1*0+0+1))

:: ---

:: +sea:  @rpb -> $up

:: 1. Read s (1 bit).
:: 2. Read RLE until different bit (n bits).
:: 3. Set e to 0.
:: 4. Get f from remaining bits and formula.

--
