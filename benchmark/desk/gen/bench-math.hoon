::  bench-math: run ONE timing cell.  +bench-math [impl door arm n]
::    impl ?(%cheb %taylor) ; door ?(%rh %rs %rd %rq) ; arm @tas (or %base) ; n @ud
::  ~>(%bout ..) inside +cell prints elapsed; returns the folded accumulator.
::  Per-call = (cell(arm) - cell(%base)) / n.  See lib/bench-cells for the guts.
::
/+  bc=bench-cells
:-  %say
|=  [* [impl=?(%cheb %taylor) door=?(%rh %rs %rd %rq) arm=@tas n=@ud ~] ~]
:-  %noun
(cell:bc impl door arm n)
