::  bench-core: shared timing loop for the numerics benchmark suite.
::
::  `time` runs a tight n-iteration loop wrapped in ~>(%bout ...), which PRINTS
::  the elapsed time ("took ms/..") as a slog and RETURNS the computed value.
::  The host driver (tools/bench_run.sh) scrapes the printed line; the returned
::  value is the folded accumulator, present only to force evaluation of every
::  call (defeats dead-code elimination).
::
::  The per-iteration work is supplied as a `step` gate |=([i=@ud acc=@] @):
::  given the loop counter i and the running accumulator, it generates a VARYING
::  input from i, calls the primitive, folds the result into acc, and returns it.
::  Because the input depends on i and the result feeds acc (which is returned),
::  neither the input nor the whole loop is a constant the runtime can memo-cache.
::
::  Per-call cost = (time(arm-step) - time(baseline-step)) / n, where the baseline
::  step does the input generation and fold but skips the primitive.
::
|%
::    +time:  run n iterations of `step`, timed by %bout; return the accumulator.
::
++  time
  |=  [n=@ud step=$-([@ud @] @)]
  ^-  @
  ~>  %bout
  =/  i=@ud   0
  =/  acc=@   `@`0
  |-  ^-  @
  ?:  =(i n)  acc
  $(i +(i), acc (step i acc))
--
