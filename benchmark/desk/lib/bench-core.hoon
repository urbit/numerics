::  bench-core: shared timing loop for the numerics benchmark suite.
::
::  `time` runs a tight loop over a PRECOMPUTED list of inputs, wrapped in
::  ~>(%bout ...), which PRINTS the elapsed time ("took ms/..") as a slog and
::  RETURNS the folded accumulator.  The host driver scrapes the printed line;
::  the returned value forces evaluation of every call (defeats dead-code
::  elimination).
::
::  CRITICAL: inputs are precomputed by the caller OUTSIDE this gate, so the
::  slow interpreted @ud->@rX conversions (sun:rd / san:rd, ~93 us/call) are
::  NOT charged to per-call cost.  Inside the timed loop the only work is the
::  arm under test plus a jetted atom-add fold -- both the input list walk
::  (O(1) head access) and the fold are cheap, so the measured time reflects
::  the arm, isolated.  Each input is a [x y] pair; y is unused (0) for the
::  single-argument arms and carries the second operand for atan2/pow/pow-n.
::
::  Per-call cost = (time(arm-list) - time(base-list)) / n, where the base list
::  uses the same inputs but the step skips the transcendental.
::
|%
::    +time:  walk the precomputed input list, timed by %bout; return the acc.
::
++  time
  |=  [xs=(list [x=@ y=@]) step=$-([[x=@ y=@] acc=@] @)]
  ^-  @
  ~>  %bout
  =/  acc=@  `@`0
  |-  ^-  @
  ?~  xs  acc
  $(xs t.xs, acc (step i.xs acc))
--
