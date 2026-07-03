::  bench-unum-grid: run the whole /lib/unum timing grid in one dojo call.
::    +bench-unum-grid n=@ud
::  For each (door, arm) it ~&-prints a [%cell door arm] label, then +cell's
::  ~>(%bout) slogs "took ..".  Per-call = (arm took - base took) / n.  Also
::  runs +fdp-cell per door (per-element = took / n).  Returns the folded
::  results to force evaluation.  Run once on the jetted build (hints on) and
::  once interpreted (hints commented) and diff the scraped times.
::
/+  uc=unum-cells
:-  %say
|=  [* [n=@ud ~] ~]
:-  %noun
::  Same transcendental coverage as math's benchmark (exp log sin cos tan atan
::  asin acos sqt cbrt pow pow-n log-2 log-10), minus atan2 (unum has no
::  atan2 arm -- see NEXT-STEPS.md), plus the arithmetic arms math's grid
::  doesn't separately test.
=/  arms=(list @tas)
  :~  %base  %add  %sub  %mul  %div  %fma  %sqt  %neg
      %exp  %log  %log-2  %log-10  %sin  %cos  %tan  %atan  %asin  %acos
      %cbrt  %pow  %pow-n  %lth
  ==
=/  doors=(list @tas)  ~[%rpb %rph %rps]
=/  rows
  %+  turn  doors
  |=  door=@tas
  %+  turn  arms
  |=  arm=@tas
  ~&  [%cell door arm]
  (cell:uc door arm n)
=/  fdps
  %+  turn  doors
  |=  door=@tas
  ~&  [%fdp door]
  (fdp-cell:uc door n)
[rows fdps]
