::  bench-unum: run ONE /lib/unum timing cell.  +bench-unum [door arm n]
::    door ?(%rpb %rph %rps) ; arm @tas (or %base) ; n @ud
::  +cell's ~>(%bout) prints elapsed; per-call = (cell(arm) - cell(%base)) / n.
::
/+  uc=unum-cells
:-  %say
|=  [* [door=?(%rpb %rph %rps) arm=@tas n=@ud ~] ~]
:-  %noun
(cell:uc door arm n)
