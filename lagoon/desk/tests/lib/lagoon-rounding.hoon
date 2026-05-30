/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-rounding -- rounding-mode behavior
::
::  Named ++test-<op>-<shape>-<kind>.  A `canon` is the reference result,
::  an `assay` is the result of the operation under test.
::
^|
|_  $:  atol=_.1e-3          :: absolute tolerance for precision of operations
        rtol=_.1e-5          :: relative tolerance for precision of operations
    ==
::  Auxiliary tools
++  is-equal
  |=  [a=ray b=ray]  ^-  tang
  ?:  =(a b)  ~
  :~  [%palm [": " ~ ~ ~] [leaf+"expected" "{<`ray`a>}"]]
      [%palm [": " ~ ~ ~] [leaf+"actual  " "{<`ray`b>}"]]
  ==
::
++  is-close
  |=  [a=ray b=ray]  ^-  tang
  ?:  (all:la (is-close:la a b [atol rtol]))  ~
  :~  [%palm [": " ~ ~ ~] [leaf+"expected" "{<a>}"]]
      [%palm [": " ~ ~ ~] [leaf+"actual  " "{<b>}"]]
  ==
::

++  test-rounding-4r-asc-z  ^-  tang
  =/  canon-cumsum-4r  `ray`[meta=[shape=~[1] bloq=4 kind=%i754 prec=~] data=0x1.457f]
  =/  assay-cumsum-4r  (cumsum:la (linspace:(lake %z) [~[11] 4 %i754 ~] [.~~1 .~~0] 11))
  %+  is-equal
    canon-cumsum-4r
  assay-cumsum-4r
::

++  test-rounding-4r-asc-d  ^-  tang
  =/  canon-cumsum-4r  `ray`[meta=[shape=~[1] bloq=4 kind=%i754 prec=~] data=0x1.457e]
  =/  assay-cumsum-4r  (cumsum:la (linspace:(lake %d) [~[11] 4 %i754 ~] [.~~1 .~~0] 11))
  %+  is-equal
    canon-cumsum-4r
  assay-cumsum-4r
::

++  test-rounding-4r-asc-u  ^-  tang
  =/  canon-cumsum-4r  `ray`[meta=[shape=~[1] bloq=4 kind=%i754 prec=~] data=0x1.457f]
  =/  assay-cumsum-4r  (cumsum:la (linspace:(lake %u) [~[11] 4 %i754 ~] [.~~1 .~~0] 11))
  %+  is-equal
    canon-cumsum-4r
  assay-cumsum-4r
::

++  test-rounding-4r-asc-n  ^-  tang
  =/  canon-cumsum-4r  `ray`[meta=[shape=~[1] bloq=4 kind=%i754 prec=~] data=0x1.457f]
  =/  assay-cumsum-4r  (cumsum:la (linspace:(lake %n) [~[11] 4 %i754 ~] [.~~1 .~~0] 11))
  %+  is-equal
    canon-cumsum-4r
  assay-cumsum-4r
::

++  test-rounding-4r-des-z  ^-  tang
  =/  canon-cumsum-4r  `ray`[meta=[shape=~[1] bloq=4 kind=%i754 prec=~] data=0x1.457d]
  =/  assay-cumsum-4r  (cumsum:la (linspace:(lake %z) [~[11] 4 %i754 ~] [.~~0 .~~1] 11))
  %+  is-equal
    canon-cumsum-4r
  assay-cumsum-4r
::

++  test-rounding-4r-des-d  ^-  tang
  =/  canon-cumsum-4r  `ray`[meta=[shape=~[1] bloq=4 kind=%i754 prec=~] data=0x1.457d]
  =/  assay-cumsum-4r  (cumsum:la (linspace:(lake %d) [~[11] 4 %i754 ~] [.~~0 .~~1] 11))
  %+  is-equal
    canon-cumsum-4r
  assay-cumsum-4r
::

++  test-rounding-4r-des-u  ^-  tang
  =/  canon-cumsum-4r  `ray`[meta=[shape=~[1] bloq=4 kind=%i754 prec=~] data=0x1.457d]
  =/  assay-cumsum-4r  (cumsum:la (linspace:(lake %u) [~[11] 4 %i754 ~] [.~~0 .~~1] 11))
  %+  is-equal
    canon-cumsum-4r
  assay-cumsum-4r
::

++  test-rounding-4r-des-n  ^-  tang
  =/  canon-cumsum-4r  `ray`[meta=[shape=~[1] bloq=4 kind=%i754 prec=~] data=0x1.457d]
  =/  assay-cumsum-4r  (cumsum:la (linspace:(lake %n) [~[11] 4 %i754 ~] [.~~0 .~~1] 11))
  %+  is-equal
    canon-cumsum-4r
  assay-cumsum-4r
::

++  test-rounding-5r-asc-z  ^-  tang
  =/  canon-cumsum-5r  `ray`[meta=[shape=~[1] bloq=5 kind=%i754 prec=~] data=0x1.40af.ffff]
  =/  assay-cumsum-5r  (cumsum:la (linspace:(lake %z) [~[11] 5 %i754 ~] [.1 .0] 11))
  %+  is-equal
    canon-cumsum-5r
  assay-cumsum-5r
::

++  test-rounding-5r-asc-d  ^-  tang
  =/  canon-cumsum-5r  `ray`[meta=[shape=~[1] bloq=5 kind=%i754 prec=~] data=0x1.40af.fffe]
  =/  assay-cumsum-5r  (cumsum:la (linspace:(lake %d) [~[11] 5 %i754 ~] [.1 .0] 11))
  %+  is-equal
    canon-cumsum-5r
  assay-cumsum-5r
::

++  test-rounding-5r-asc-u  ^-  tang
  =/  canon-cumsum-5r  `ray`[meta=[shape=~[1] bloq=5 kind=%i754 prec=~] data=0x1.40b0.0000]
  =/  assay-cumsum-5r  (cumsum:la (linspace:(lake %u) [~[11] 5 %i754 ~] [.1 .0] 11))
  %+  is-equal
    canon-cumsum-5r
  assay-cumsum-5r
::

++  test-rounding-5r-asc-n  ^-  tang
  =/  canon-cumsum-5r  `ray`[meta=[shape=~[1] bloq=5 kind=%i754 prec=~] data=0x1.40af.ffff]
  =/  assay-cumsum-5r  (cumsum:la (linspace:(lake %n) [~[11] 5 %i754 ~] [.1 .0] 11))
  %+  is-equal
    canon-cumsum-5r
  assay-cumsum-5r
::

++  test-rounding-5r-des-z  ^-  tang
  =/  canon-cumsum-5r  `ray`[meta=[shape=~[1] bloq=5 kind=%i754 prec=~] data=0x1.40af.fffc]
  =/  assay-cumsum-5r  (cumsum:la (linspace:(lake %z) [~[11] 5 %i754 ~] [.0 .1] 11))
  %+  is-equal
    canon-cumsum-5r
  assay-cumsum-5r
::

++  test-rounding-5r-des-d  ^-  tang
  =/  canon-cumsum-5r  `ray`[meta=[shape=~[1] bloq=5 kind=%i754 prec=~] data=0x1.40af.fffc]
  =/  assay-cumsum-5r  (cumsum:la (linspace:(lake %d) [~[11] 5 %i754 ~] [.0 .1] 11))
  %+  is-equal
    canon-cumsum-5r
  assay-cumsum-5r
::

++  test-rounding-5r-des-u  ^-  tang
  =/  canon-cumsum-5r  `ray`[meta=[shape=~[1] bloq=5 kind=%i754 prec=~] data=0x1.40af.fffd]
  =/  assay-cumsum-5r  (cumsum:la (linspace:(lake %u) [~[11] 5 %i754 ~] [.0 .1] 11))
  %+  is-equal
    canon-cumsum-5r
  assay-cumsum-5r
::

++  test-rounding-5r-des-n  ^-  tang
  =/  canon-cumsum-5r  `ray`[meta=[shape=~[1] bloq=5 kind=%i754 prec=~] data=0x1.40af.fffd]
  =/  assay-cumsum-5r  (cumsum:la (linspace:(lake %n) [~[11] 5 %i754 ~] [.0 .1] 11))
  %+  is-equal
    canon-cumsum-5r
  assay-cumsum-5r
::

++  test-rounding-6r-asc-z  ^-  tang
  =/  canon-cumsum-6r  `ray`[meta=[shape=~[1] bloq=6 kind=%i754 prec=~] data=0x1.4015.ffff.ffff.ffff]
  =/  assay-cumsum-6r  (cumsum:la (linspace:(lake %z) [~[11] 6 %i754 ~] [.~1 .~0] 11))
  %+  is-equal
    canon-cumsum-6r
  assay-cumsum-6r
::

++  test-rounding-6r-asc-d  ^-  tang
  =/  canon-cumsum-6r  `ray`[meta=[shape=~[1] bloq=6 kind=%i754 prec=~] data=0x1.4015.ffff.ffff.fffe]
  =/  assay-cumsum-6r  (cumsum:la (linspace:(lake %d) [~[11] 6 %i754 ~] [.~1 .~0] 11))
  %+  is-equal
    canon-cumsum-6r
  assay-cumsum-6r
::

++  test-rounding-6r-asc-u  ^-  tang
  =/  canon-cumsum-6r  `ray`[meta=[shape=~[1] bloq=6 kind=%i754 prec=~] data=0x1.4016.0000.0000.0000]
  =/  assay-cumsum-6r  (cumsum:la (linspace:(lake %u) [~[11] 6 %i754 ~] [.~1 .~0] 11))
  %+  is-equal
    canon-cumsum-6r
  assay-cumsum-6r
::

++  test-rounding-6r-asc-n  ^-  tang
  =/  canon-cumsum-6r  `ray`[meta=[shape=~[1] bloq=6 kind=%i754 prec=~] data=0x1.4015.ffff.ffff.ffff]
  =/  assay-cumsum-6r  (cumsum:la (linspace:(lake %n) [~[11] 6 %i754 ~] [.~1 .~0] 11))
  %+  is-equal
    canon-cumsum-6r
  assay-cumsum-6r
::

++  test-rounding-6r-des-z  ^-  tang
  =/  canon-cumsum-6r  `ray`[meta=[shape=~[1] bloq=6 kind=%i754 prec=~] data=0x1.4015.ffff.ffff.fffc]
  =/  assay-cumsum-6r  (cumsum:la (linspace:(lake %z) [~[11] 6 %i754 ~] [.~0 .~1] 11))
  %+  is-equal
    canon-cumsum-6r
  assay-cumsum-6r
::

++  test-rounding-6r-des-d  ^-  tang
  =/  canon-cumsum-6r  `ray`[meta=[shape=~[1] bloq=6 kind=%i754 prec=~] data=0x1.4015.ffff.ffff.fffc]
  =/  assay-cumsum-6r  (cumsum:la (linspace:(lake %d) [~[11] 6 %i754 ~] [.~0 .~1] 11))
  %+  is-equal
    canon-cumsum-6r
  assay-cumsum-6r
::

++  test-rounding-6r-des-u  ^-  tang
  =/  canon-cumsum-6r  `ray`[meta=[shape=~[1] bloq=6 kind=%i754 prec=~] data=0x1.4015.ffff.ffff.fffe]
  =/  assay-cumsum-6r  (cumsum:la (linspace:(lake %u) [~[11] 6 %i754 ~] [.~0 .~1] 11))
  %+  is-equal
    canon-cumsum-6r
  assay-cumsum-6r
::

++  test-rounding-6r-des-n  ^-  tang
  =/  canon-cumsum-6r  `ray`[meta=[shape=~[1] bloq=6 kind=%i754 prec=~] data=0x1.4015.ffff.ffff.fffd]
  =/  assay-cumsum-6r  (cumsum:la (linspace:(lake %n) [~[11] 6 %i754 ~] [.~0 .~1] 11))
  %+  is-equal
    canon-cumsum-6r
  assay-cumsum-6r
::

++  test-rounding-7r-asc-z  ^-  tang
  =/  canon-cumsum-7r  `ray`[meta=[shape=~[1] bloq=7 kind=%i754 prec=~] data=0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.ffff]
  =/  assay-cumsum-7r  (cumsum:la (linspace:(lake %z) [~[11] 7 %i754 ~] [.~~~1 .~~~0] 11))
  %+  is-equal
    canon-cumsum-7r
  assay-cumsum-7r
::

++  test-rounding-7r-asc-d  ^-  tang
  =/  canon-cumsum-7r  `ray`[meta=[shape=~[1] bloq=7 kind=%i754 prec=~] data=0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.fffe]
  =/  assay-cumsum-7r  (cumsum:la (linspace:(lake %d) [~[11] 7 %i754 ~] [.~~~1 .~~~0] 11))
  %+  is-equal
    canon-cumsum-7r
  assay-cumsum-7r
::

++  test-rounding-7r-asc-u  ^-  tang
  =/  canon-cumsum-7r  `ray`[meta=[shape=~[1] bloq=7 kind=%i754 prec=~] data=0x1.4001.6000.0000.0000.0000.0000.0000.0000]
  =/  assay-cumsum-7r  (cumsum:la (linspace:(lake %u) [~[11] 7 %i754 ~] [.~~~1 .~~~0] 11))
  %+  is-equal
    canon-cumsum-7r
  assay-cumsum-7r
::

++  test-rounding-7r-asc-n  ^-  tang
  =/  canon-cumsum-7r  `ray`[meta=[shape=~[1] bloq=7 kind=%i754 prec=~] data=0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.ffff]
  =/  assay-cumsum-7r  (cumsum:la (linspace:(lake %n) [~[11] 7 %i754 ~] [.~~~1 .~~~0] 11))
  %+  is-equal
    canon-cumsum-7r
  assay-cumsum-7r
::

++  test-rounding-7r-des-z  ^-  tang
  =/  canon-cumsum-7r  `ray`[meta=[shape=~[1] bloq=7 kind=%i754 prec=~] data=0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.fffc]
  =/  assay-cumsum-7r  (cumsum:la (linspace:(lake %z) [~[11] 7 %i754 ~] [.~~~0 .~~~1] 11))
  %+  is-equal
    canon-cumsum-7r
  assay-cumsum-7r
::

++  test-rounding-7r-des-d  ^-  tang
  =/  canon-cumsum-7r  `ray`[meta=[shape=~[1] bloq=7 kind=%i754 prec=~] data=0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.fffc]
  =/  assay-cumsum-7r  (cumsum:la (linspace:(lake %d) [~[11] 7 %i754 ~] [.~~~0 .~~~1] 11))
  %+  is-equal
    canon-cumsum-7r
  assay-cumsum-7r
::

++  test-rounding-7r-des-u  ^-  tang
  =/  canon-cumsum-7r  `ray`[meta=[shape=~[1] bloq=7 kind=%i754 prec=~] data=0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.fffe]
  =/  assay-cumsum-7r  (cumsum:la (linspace:(lake %u) [~[11] 7 %i754 ~] [.~~~0 .~~~1] 11))
  %+  is-equal
    canon-cumsum-7r
  assay-cumsum-7r
::

++  test-rounding-7r-des-n  ^-  tang
  =/  canon-cumsum-7r  `ray`[meta=[shape=~[1] bloq=7 kind=%i754 prec=~] data=0x1.4001.5fff.ffff.ffff.ffff.ffff.ffff.fffd]
  =/  assay-cumsum-7r  (cumsum:la (linspace:(lake %n) [~[11] 7 %i754 ~] [.~~~0 .~~~1] 11))
  %+  is-equal
    canon-cumsum-7r
  assay-cumsum-7r
--
