/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-rounding -- rounding-mode behavior
::
::  cumsum of an 11-point linspace, built under each rounding mode, must
::  match a mode-specific reference value.  Table-driven: one row per
::  (precision, direction, mode), one iterator arm, instead of 32
::  near-identical hand-written arms.
::
^|
|_  $:  atol=_.1e-3          :: absolute tolerance for precision of operations
        rtol=_.1e-5          :: relative tolerance for precision of operations
    ==
++  is-equal
  |=  [a=ray b=ray]  ^-  tang
  ?:  =(a b)  ~
  :~  [%palm [": " ~ ~ ~] [leaf+"expected" "{<`ray`a>}"]]
      [%palm [": " ~ ~ ~] [leaf+"actual  " "{<`ray`b>}"]]
  ==
::
::  [bloq mode lo hi canon] -- linspace runs [lo hi] over 11 points under
::  `mode`, then cumsum (always default mode); canon is the summed scalar.
++  rnd-table
  ^-  (list [bloq=@ mode=?(%n %u %d %z) lo=@ hi=@ canon=@])
  :~  [4 %z .~~1 .~~0 0x1.457f]
      [4 %d .~~1 .~~0 0x1.457e]
      [4 %u .~~1 .~~0 0x1.457f]
      [4 %n .~~1 .~~0 0x1.457f]
      [4 %z .~~0 .~~1 0x1.457d]
      [4 %d .~~0 .~~1 0x1.457d]
      [4 %u .~~0 .~~1 0x1.457d]
      [4 %n .~~0 .~~1 0x1.457d]
      [5 %z .1 .0 0x1.40af.ffff]
      [5 %d .1 .0 0x1.40af.fffe]
      [5 %u .1 .0 0x1.40b0.0000]
      [5 %n .1 .0 0x1.40af.ffff]
      [5 %z .0 .1 0x1.40af.fffc]
      [5 %d .0 .1 0x1.40af.fffc]
      [5 %u .0 .1 0x1.40af.fffd]
      [5 %n .0 .1 0x1.40af.fffd]
      [6 %z .~1 .~0 0x1.4015.ffff.ffff.ffff]
      [6 %d .~1 .~0 0x1.4015.ffff.ffff.fffe]
      [6 %u .~1 .~0 0x1.4016.0000.0000.0000]
      [6 %n .~1 .~0 0x1.4015.ffff.ffff.ffff]
      [6 %z .~0 .~1 0x1.4015.ffff.ffff.fffc]
      [6 %d .~0 .~1 0x1.4015.ffff.ffff.fffc]
      [6 %u .~0 .~1 0x1.4015.ffff.ffff.fffe]
      [6 %n .~0 .~1 0x1.4015.ffff.ffff.fffd]
      [7 %z .~~~1 .~~~0 0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.ffff]
      [7 %d .~~~1 .~~~0 0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.fffe]
      [7 %u .~~~1 .~~~0 0x1.4001.6000.0000.0000.0000.0000.0000.0000]
      [7 %n .~~~1 .~~~0 0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.ffff]
      [7 %z .~~~0 .~~~1 0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.fffc]
      [7 %d .~~~0 .~~~1 0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.fffc]
      [7 %u .~~~0 .~~~1 0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.fffe]
      [7 %n .~~~0 .~~~1 0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.fffd]
  ==
++  test-rounding  ^-  tang
  %-  zing
  %+  turn  rnd-table
  |=  [bloq=@ mode=?(%n %u %d %z) lo=@ hi=@ canon=@]
  ^-  tang
  %+  is-equal
    `ray`[[~[1] bloq %i754 ~] canon]
  (cumsum:la (linspace:(lake mode) [~[11] bloq %i754 ~] [lo hi] 11))
--
