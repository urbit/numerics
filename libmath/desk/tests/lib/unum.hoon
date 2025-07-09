  ::  /tests/lib/unum
::::
::    Posits
::
/+  math,
    *test,
    unum
^|
|%
::
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
      !>  `@rpb`0b1111.1111
      !>  neg-huge:rpb:unum
    ::  tiny
    %+  expect-eq
      !>  `@rpb`0b1
      !>  tiny:rpb:unum
  ==
::
++  test-sea-rpb  ^-  tang
  ;:  weld
    ::  0
    %+  expect-eq
      !>  [%z 3 ~]
      !>  (sea:rpb:unum zero:rpb:unum)
    ::  1
    %+  expect-eq
      !>  [%p 3 s=%.y r=--0 e=--0 f=0]
      !>  (sea:rpb:unum one:rpb:unum)
  ==
::
:: ++  test-abs  ^-  tang
::   ;:  weld
::     %+  expect-near
::       .1
::       (abs:rs:saloon .-1)
::     %+  expect-near
::       .1
::       (abs:rs:saloon .1)
::     %+  expect-near
::       .120
::       (abs:rs:saloon .-120)
::     ==  
--
