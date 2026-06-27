::  bench-grid: run the WHOLE timing grid in ONE dojo invocation.
::    +bench-grid [impl sets nrd nrs nrh nrq]      impl ?(%cheb %taylor)
::  Per-door iteration counts (nrd/nrs/nrh/nrq) so the jetted grid can use 100k
::  everywhere while the INTERPRETED grid winds @rh/@rq way down -- their base
::  kernel ops aren't jetted on this binary, so interpreted f16/f128 is ~1000x
::  slower and would otherwise hang the single-command grid.
::  For every (door, arm, set) it slogs a marker `[%cell door arm set n]` then runs
::  the cell, whose ~>(%bout) slogs `took ..`.  The host scrapes marker+took pairs
::  (n is read from the marker).  Each cell is mule-guarded: a crash slogs
::  `[%fail door arm]` and the grid continues.  Ends with `[%grid-done]`.  Returns ~.
::
/+  bc=bench-cells
:-  %say
|=  [* [impl=?(%cheb %taylor) sets=@ud nrd=@ud nrs=@ud nrh=@ud nrq=@ud scope=?(%all %two) ~] ~]
:-  %noun
::  scope %two = @rd/@rs only (interpreted @rh/@rq are infeasible: un-jetted base ops).
=/  doors=(list @tas)  ?:(=(scope %two) ~[%rd %rs] ~[%rd %rs %rh %rq])
=/  arms=(list @tas)
  :~  %base  %exp  %log  %sin  %cos  %tan  %atan  %atan2
      %asin  %acos  %sqt  %cbt  %pow  %pow-n  %log-2  %log-10
  ==
=/  triples
  ^-  (list [@tas @tas @ud])
  %-  zing
  %+  turn  doors
  |=  d=@tas  ^-  (list [@tas @tas @ud])
  %-  zing
  %+  turn  arms
  |=  a=@tas  ^-  (list [@tas @tas @ud])
  %+  turn  (gulf 0 (dec sets))
  |=  s=@ud  ^-  [@tas @tas @ud]
  [d a s]
|-  ^-  ~
?~  triples  ~&(> [%grid-done] ~)
=*  d  -.i.triples
=*  a  +<.i.triples
=*  s  +>.i.triples
=/  cn=@ud
  ?+  d  nrd
    %rs  nrs
    %rh  nrh
    %rq  nrq
  ==
~&  >  [%cell d a s cn]
=/  res  (mule |.(`@`(cell:bc impl d a cn)))
?:  ?=(%& -.res)
  $(triples t.triples)
~&(>>> [%fail d a] $(triples t.triples))
