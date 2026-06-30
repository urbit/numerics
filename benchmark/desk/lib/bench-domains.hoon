::  bench-domains: per-arm input domains + counter->input generators.
::
::  Domains are stored as INTEGER bounds [lo hi denom] meaning the real interval
::  [lo/denom, hi/denom], so one table serves all four float widths; the float is
::  produced per-door from the base door's `san`/`sun` (jetted) ops.  `span` is the
::  number of distinct samples to cycle through (capped for @rh, whose ULP is coarse).
::
::  `kind`:  %one  single-arg arm;  %two  second operand stepped too (atan2);
::          %pwr  second operand is a fixed float exponent (pow);
::          %pwn  second operand is a small integer exponent cycled 2..8 (pow-n).
::
::  Domains are DOOR-AWARE: @rh (f16) has a tiny dynamic range (~6e-5 .. 65504),
::  so the f64-scale intervals push exp/pow-n into under/overflow, where the
::  INTERPRETED transcendentals churn (>300x slower).  For @rh those arms get
::  tighter, in-range intervals; @rs/@rd/@rq use the base table (f32+ range is
::  ample).  Accuracy/speed are still measured over each width's sensible domain.
::
|%
+$  dom  [lo=@sd hi=@sd den=@ud span=@ud kind=?(%one %two %pwr %pwn)]
::    +get:  domain for (arm, impl, door).  Taylor uses tighter "safe" intervals.
::
++  get
  |=  [arm=@tas impl=?(%cheb %taylor) door=@tas]
  ^-  dom
  =/  taylor  =(impl %taylor)
  =/  base=dom
    ?+  arm  ~|([%bad-arm arm] !!)
      %exp     ?:(taylor [-5 --5 1 10.000 %one] [-20 --20 1 40.000 %one])
      %log     ?:(taylor [--1 --100 10 10.000 %one] [--1 --2.000 10 20.000 %one])
      %sin     ?:(taylor [-31 --31 10 6.283 %one] [-100 --100 1 40.000 %one])
      %cos     ?:(taylor [-31 --31 10 6.283 %one] [-100 --100 1 40.000 %one])
      %tan     ?:(taylor [-1 --1 1 4.000 %one] [-3 --3 2 6.000 %one])
      %atan    ?:(taylor [-5 --5 1 10.000 %one] [-50 --50 1 40.000 %one])
      %atan2   [-10 --10 1 20.000 %two]
      %asin    ?:(taylor [-99 --99 100 4.000 %one] [-1 --1 1 4.000 %one])
      %acos    ?:(taylor [-99 --99 100 4.000 %one] [-1 --1 1 4.000 %one])
      %sqt     ?:(taylor [--1 --10.000 100 20.000 %one] [--1 --20.000 100 20.000 %one])
      %cbt     ?:(taylor [--1 --10.000 100 20.000 %one] [-20.000 --20.000 100 40.000 %one])
      %pow     [--1 --10 1 10.000 %pwr]
      %pow-n   [--1 --10 1 10.000 %pwn]
      %log-2   ?:(taylor [--1 --100 10 10.000 %one] [--1 --2.000 10 20.000 %one])
      %log-10  ?:(taylor [--1 --100 10 10.000 %one] [--1 --2.000 10 20.000 %one])
    ==
  ::  @rh f16-safe overrides for the arms whose f64-scale range under/overflows.
  ?.  =(door %rh)  base
  ?+  arm  base
    %exp     [-8 --8 1 8.000 %one]        :: exp(+-8) in [3.4e-4, 2981], in-range
    %sin     [-12 --12 1 24.000 %one]     :: ~2 periods; avoids huge range-reduction
    %cos     [-12 --12 1 24.000 %one]
    %atan    [-16 --16 1 32.000 %one]
    %pow-n   [--1 --3 1 3.000 %pwn]       :: 3^8 = 6561 < 65504, no overflow
    %log     [--1 --400 10 20.000 %one]   :: (0.1, 40]; log stays small either way
    %log-2   [--1 --400 10 20.000 %one]
    %log-10  [--1 --400 10 20.000 %one]
  ==
--
