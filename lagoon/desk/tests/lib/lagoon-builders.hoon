/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-builders -- array builders
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
++  test-eye-1x1-4r  ^-  tang
  =/  input-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0]]])
  =/  canon-eye-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0]]])
  =/  assay-eye-1x1-4r  (eye:la meta.input-ones-1x1-4r)
  %+  is-equal
    canon-eye-1x1-4r
  assay-eye-1x1-4r
::

++  test-eye-2x2-4r  ^-  tang
  =/  input-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-eye-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~0.0] ~[.~~0.0 .~~1.0]]])
  =/  assay-eye-2x2-4r  (eye:la meta.input-ones-2x2-4r)
  %+  is-equal
    canon-eye-2x2-4r
  assay-eye-2x2-4r
::

++  test-eye-3x3-4r  ^-  tang
  =/  input-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-eye-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~0.0 .~~0.0] ~[.~~0.0 .~~1.0 .~~0.0] ~[.~~0.0 .~~0.0 .~~1.0]]])
  =/  assay-eye-3x3-4r  (eye:la meta.input-ones-3x3-4r)
  %+  is-equal
    canon-eye-3x3-4r
  assay-eye-3x3-4r
::

++  test-eye-1x1-5r  ^-  tang
  =/  input-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0]]])
  =/  canon-eye-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0]]])
  =/  assay-eye-1x1-5r  (eye:la meta.input-ones-1x1-5r)
  %+  is-equal
    canon-eye-1x1-5r
  assay-eye-1x1-5r
::

++  test-eye-2x2-5r  ^-  tang
  =/  input-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-eye-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .0.0] ~[.0.0 .1.0]]])
  =/  assay-eye-2x2-5r  (eye:la meta.input-ones-2x2-5r)
  %+  is-equal
    canon-eye-2x2-5r
  assay-eye-2x2-5r
::

++  test-eye-3x3-5r  ^-  tang
  =/  input-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-eye-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .0.0 .0.0] ~[.0.0 .1.0 .0.0] ~[.0.0 .0.0 .1.0]]])
  =/  assay-eye-3x3-5r  (eye:la meta.input-ones-3x3-5r)
  %+  is-equal
    canon-eye-3x3-5r
  assay-eye-3x3-5r
::

++  test-eye-1x1-6r  ^-  tang
  =/  input-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0]]])
  =/  canon-eye-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0]]])
  =/  assay-eye-1x1-6r  (eye:la meta.input-ones-1x1-6r)
  %+  is-equal
    canon-eye-1x1-6r
  assay-eye-1x1-6r
::

++  test-eye-2x2-6r  ^-  tang
  =/  input-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-eye-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~0.0] ~[.~0.0 .~1.0]]])
  =/  assay-eye-2x2-6r  (eye:la meta.input-ones-2x2-6r)
  %+  is-equal
    canon-eye-2x2-6r
  assay-eye-2x2-6r
::

++  test-eye-3x3-6r  ^-  tang
  =/  input-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-eye-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~0.0 .~0.0] ~[.~0.0 .~1.0 .~0.0] ~[.~0.0 .~0.0 .~1.0]]])
  =/  assay-eye-3x3-6r  (eye:la meta.input-ones-3x3-6r)
  %+  is-equal
    canon-eye-3x3-6r
  assay-eye-3x3-6r
::

++  test-eye-1x1-7r  ^-  tang
  =/  input-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0]]])
  =/  canon-eye-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0]]])
  =/  assay-eye-1x1-7r  (eye:la meta.input-ones-1x1-7r)
  %+  is-equal
    canon-eye-1x1-7r
  assay-eye-1x1-7r
::

++  test-eye-2x2-7r  ^-  tang
  =/  input-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-eye-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~0.0] ~[.~~~0.0 .~~~1.0]]])
  =/  assay-eye-2x2-7r  (eye:la meta.input-ones-2x2-7r)
  %+  is-equal
    canon-eye-2x2-7r
  assay-eye-2x2-7r
::

++  test-eye-3x3-7r  ^-  tang
  =/  input-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-eye-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~1.0 .~~~0.0] ~[.~~~0.0 .~~~0.0 .~~~1.0]]])
  =/  assay-eye-3x3-7r  (eye:la meta.input-ones-3x3-7r)
  %+  is-equal
    canon-eye-3x3-7r
  assay-eye-3x3-7r
::

++  test-eye-1x1-3u  ^-  tang
  =/  input-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-eye-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[1]]])
  =/  assay-eye-1x1-3u  (eye:la meta.input-ones-1x1-3u)
  %+  is-equal
    canon-eye-1x1-3u
  assay-eye-1x1-3u
::

++  test-eye-2x2-3u  ^-  tang
  =/  input-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-eye-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint prec=~] baum=~[~[1 0] ~[0 1]]])
  =/  assay-eye-2x2-3u  (eye:la meta.input-ones-2x2-3u)
  %+  is-equal
    canon-eye-2x2-3u
  assay-eye-2x2-3u
::

++  test-eye-3x3-3u  ^-  tang
  =/  input-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-eye-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint prec=~] baum=~[~[1 0 0] ~[0 1 0] ~[0 0 1]]])
  =/  assay-eye-3x3-3u  (eye:la meta.input-ones-3x3-3u)
  %+  is-equal
    canon-eye-3x3-3u
  assay-eye-3x3-3u
::

++  test-eye-1x1-4u  ^-  tang
  =/  input-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-eye-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[1]]])
  =/  assay-eye-1x1-4u  (eye:la meta.input-ones-1x1-4u)
  %+  is-equal
    canon-eye-1x1-4u
  assay-eye-1x1-4u
::

++  test-eye-2x2-4u  ^-  tang
  =/  input-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-eye-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint prec=~] baum=~[~[1 0] ~[0 1]]])
  =/  assay-eye-2x2-4u  (eye:la meta.input-ones-2x2-4u)
  %+  is-equal
    canon-eye-2x2-4u
  assay-eye-2x2-4u
::

++  test-eye-3x3-4u  ^-  tang
  =/  input-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-eye-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint prec=~] baum=~[~[1 0 0] ~[0 1 0] ~[0 0 1]]])
  =/  assay-eye-3x3-4u  (eye:la meta.input-ones-3x3-4u)
  %+  is-equal
    canon-eye-3x3-4u
  assay-eye-3x3-4u
::

++  test-eye-1x1-5u  ^-  tang
  =/  input-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-eye-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[1]]])
  =/  assay-eye-1x1-5u  (eye:la meta.input-ones-1x1-5u)
  %+  is-equal
    canon-eye-1x1-5u
  assay-eye-1x1-5u
::

++  test-eye-2x2-5u  ^-  tang
  =/  input-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-eye-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint prec=~] baum=~[~[1 0] ~[0 1]]])
  =/  assay-eye-2x2-5u  (eye:la meta.input-ones-2x2-5u)
  %+  is-equal
    canon-eye-2x2-5u
  assay-eye-2x2-5u
::

++  test-eye-3x3-5u  ^-  tang
  =/  input-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-eye-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint prec=~] baum=~[~[1 0 0] ~[0 1 0] ~[0 0 1]]])
  =/  assay-eye-3x3-5u  (eye:la meta.input-ones-3x3-5u)
  %+  is-equal
    canon-eye-3x3-5u
  assay-eye-3x3-5u
::

++  test-eye-1x1-6u  ^-  tang
  =/  input-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-eye-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[1]]])
  =/  assay-eye-1x1-6u  (eye:la meta.input-ones-1x1-6u)
  %+  is-equal
    canon-eye-1x1-6u
  assay-eye-1x1-6u
::

++  test-eye-2x2-6u  ^-  tang
  =/  input-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-eye-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint prec=~] baum=~[~[1 0] ~[0 1]]])
  =/  assay-eye-2x2-6u  (eye:la meta.input-ones-2x2-6u)
  %+  is-equal
    canon-eye-2x2-6u
  assay-eye-2x2-6u
::

++  test-eye-3x3-6u  ^-  tang
  =/  input-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-eye-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint prec=~] baum=~[~[1 0 0] ~[0 1 0] ~[0 0 1]]])
  =/  assay-eye-3x3-6u  (eye:la meta.input-ones-3x3-6u)
  %+  is-equal
    canon-eye-3x3-6u
  assay-eye-3x3-6u
::

++  test-zeros-1x1-4r  ^-  tang
  =/  input-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0]]])
  =/  canon-zeros-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0]]])
  =/  assay-zeros-1x1-4r  (zeros:la meta.input-ones-1x1-4r)
  %+  is-equal
    canon-zeros-1x1-4r
  assay-zeros-1x1-4r
::

++  test-zeros-1x2-4r  ^-  tang
  =/  input-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  canon-zeros-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0 .~~0.0]]])
  =/  assay-zeros-1x2-4r  (zeros:la meta.input-ones-1x2-4r)
  %+  is-equal
    canon-zeros-1x2-4r
  assay-zeros-1x2-4r
::

++  test-zeros-1x3-4r  ^-  tang
  =/  input-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-zeros-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0 .~~0.0 .~~0.0]]])
  =/  assay-zeros-1x3-4r  (zeros:la meta.input-ones-1x3-4r)
  %+  is-equal
    canon-zeros-1x3-4r
  assay-zeros-1x3-4r
::

++  test-zeros-2x1-4r  ^-  tang
  =/  input-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  canon-zeros-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0] ~[.~~0.0]]])
  =/  assay-zeros-2x1-4r  (zeros:la meta.input-ones-2x1-4r)
  %+  is-equal
    canon-zeros-2x1-4r
  assay-zeros-2x1-4r
::

++  test-zeros-2x2-4r  ^-  tang
  =/  input-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-zeros-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0 .~~0.0] ~[.~~0.0 .~~0.0]]])
  =/  assay-zeros-2x2-4r  (zeros:la meta.input-ones-2x2-4r)
  %+  is-equal
    canon-zeros-2x2-4r
  assay-zeros-2x2-4r
::

++  test-zeros-2x3-4r  ^-  tang
  =/  input-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-zeros-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0 .~~0.0 .~~0.0] ~[.~~0.0 .~~0.0 .~~0.0]]])
  =/  assay-zeros-2x3-4r  (zeros:la meta.input-ones-2x3-4r)
  %+  is-equal
    canon-zeros-2x3-4r
  assay-zeros-2x3-4r
::

++  test-zeros-3x1-4r  ^-  tang
  =/  input-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  canon-zeros-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0] ~[.~~0.0] ~[.~~0.0]]])
  =/  assay-zeros-3x1-4r  (zeros:la meta.input-ones-3x1-4r)
  %+  is-equal
    canon-zeros-3x1-4r
  assay-zeros-3x1-4r
::

++  test-zeros-3x2-4r  ^-  tang
  =/  input-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-zeros-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0 .~~0.0] ~[.~~0.0 .~~0.0] ~[.~~0.0 .~~0.0]]])
  =/  assay-zeros-3x2-4r  (zeros:la meta.input-ones-3x2-4r)
  %+  is-equal
    canon-zeros-3x2-4r
  assay-zeros-3x2-4r
::

++  test-zeros-3x3-4r  ^-  tang
  =/  input-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-zeros-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0 .~~0.0 .~~0.0] ~[.~~0.0 .~~0.0 .~~0.0] ~[.~~0.0 .~~0.0 .~~0.0]]])
  =/  assay-zeros-3x3-4r  (zeros:la meta.input-ones-3x3-4r)
  %+  is-equal
    canon-zeros-3x3-4r
  assay-zeros-3x3-4r
::

++  test-zeros-1x1-5r  ^-  tang
  =/  input-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0]]])
  =/  canon-zeros-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0]]])
  =/  assay-zeros-1x1-5r  (zeros:la meta.input-ones-1x1-5r)
  %+  is-equal
    canon-zeros-1x1-5r
  assay-zeros-1x1-5r
::

++  test-zeros-1x2-5r  ^-  tang
  =/  input-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0]]])
  =/  canon-zeros-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0 .0.0]]])
  =/  assay-zeros-1x2-5r  (zeros:la meta.input-ones-1x2-5r)
  %+  is-equal
    canon-zeros-1x2-5r
  assay-zeros-1x2-5r
::

++  test-zeros-1x3-5r  ^-  tang
  =/  input-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  canon-zeros-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0 .0.0 .0.0]]])
  =/  assay-zeros-1x3-5r  (zeros:la meta.input-ones-1x3-5r)
  %+  is-equal
    canon-zeros-1x3-5r
  assay-zeros-1x3-5r
::

++  test-zeros-2x1-5r  ^-  tang
  =/  input-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  canon-zeros-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0] ~[.0.0]]])
  =/  assay-zeros-2x1-5r  (zeros:la meta.input-ones-2x1-5r)
  %+  is-equal
    canon-zeros-2x1-5r
  assay-zeros-2x1-5r
::

++  test-zeros-2x2-5r  ^-  tang
  =/  input-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-zeros-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0 .0.0] ~[.0.0 .0.0]]])
  =/  assay-zeros-2x2-5r  (zeros:la meta.input-ones-2x2-5r)
  %+  is-equal
    canon-zeros-2x2-5r
  assay-zeros-2x2-5r
::

++  test-zeros-2x3-5r  ^-  tang
  =/  input-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-zeros-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0 .0.0 .0.0] ~[.0.0 .0.0 .0.0]]])
  =/  assay-zeros-2x3-5r  (zeros:la meta.input-ones-2x3-5r)
  %+  is-equal
    canon-zeros-2x3-5r
  assay-zeros-2x3-5r
::

++  test-zeros-3x1-5r  ^-  tang
  =/  input-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  canon-zeros-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0] ~[.0.0] ~[.0.0]]])
  =/  assay-zeros-3x1-5r  (zeros:la meta.input-ones-3x1-5r)
  %+  is-equal
    canon-zeros-3x1-5r
  assay-zeros-3x1-5r
::

++  test-zeros-3x2-5r  ^-  tang
  =/  input-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-zeros-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0 .0.0] ~[.0.0 .0.0] ~[.0.0 .0.0]]])
  =/  assay-zeros-3x2-5r  (zeros:la meta.input-ones-3x2-5r)
  %+  is-equal
    canon-zeros-3x2-5r
  assay-zeros-3x2-5r
::

++  test-zeros-3x3-5r  ^-  tang
  =/  input-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-zeros-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0 .0.0 .0.0] ~[.0.0 .0.0 .0.0] ~[.0.0 .0.0 .0.0]]])
  =/  assay-zeros-3x3-5r  (zeros:la meta.input-ones-3x3-5r)
  %+  is-equal
    canon-zeros-3x3-5r
  assay-zeros-3x3-5r
::

++  test-zeros-1x1-6r  ^-  tang
  =/  input-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0]]])
  =/  canon-zeros-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0]]])
  =/  assay-zeros-1x1-6r  (zeros:la meta.input-ones-1x1-6r)
  %+  is-equal
    canon-zeros-1x1-6r
  assay-zeros-1x1-6r
::

++  test-zeros-1x2-6r  ^-  tang
  =/  input-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0]]])
  =/  canon-zeros-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0 .~0.0]]])
  =/  assay-zeros-1x2-6r  (zeros:la meta.input-ones-1x2-6r)
  %+  is-equal
    canon-zeros-1x2-6r
  assay-zeros-1x2-6r
::

++  test-zeros-1x3-6r  ^-  tang
  =/  input-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-zeros-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0 .~0.0 .~0.0]]])
  =/  assay-zeros-1x3-6r  (zeros:la meta.input-ones-1x3-6r)
  %+  is-equal
    canon-zeros-1x3-6r
  assay-zeros-1x3-6r
::

++  test-zeros-2x1-6r  ^-  tang
  =/  input-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  canon-zeros-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0] ~[.~0.0]]])
  =/  assay-zeros-2x1-6r  (zeros:la meta.input-ones-2x1-6r)
  %+  is-equal
    canon-zeros-2x1-6r
  assay-zeros-2x1-6r
::

++  test-zeros-2x2-6r  ^-  tang
  =/  input-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-zeros-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0 .~0.0] ~[.~0.0 .~0.0]]])
  =/  assay-zeros-2x2-6r  (zeros:la meta.input-ones-2x2-6r)
  %+  is-equal
    canon-zeros-2x2-6r
  assay-zeros-2x2-6r
::

++  test-zeros-2x3-6r  ^-  tang
  =/  input-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-zeros-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0 .~0.0 .~0.0] ~[.~0.0 .~0.0 .~0.0]]])
  =/  assay-zeros-2x3-6r  (zeros:la meta.input-ones-2x3-6r)
  %+  is-equal
    canon-zeros-2x3-6r
  assay-zeros-2x3-6r
::

++  test-zeros-3x1-6r  ^-  tang
  =/  input-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  canon-zeros-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0] ~[.~0.0] ~[.~0.0]]])
  =/  assay-zeros-3x1-6r  (zeros:la meta.input-ones-3x1-6r)
  %+  is-equal
    canon-zeros-3x1-6r
  assay-zeros-3x1-6r
::

++  test-zeros-3x2-6r  ^-  tang
  =/  input-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-zeros-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0 .~0.0] ~[.~0.0 .~0.0] ~[.~0.0 .~0.0]]])
  =/  assay-zeros-3x2-6r  (zeros:la meta.input-ones-3x2-6r)
  %+  is-equal
    canon-zeros-3x2-6r
  assay-zeros-3x2-6r
::

++  test-zeros-3x3-6r  ^-  tang
  =/  input-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-zeros-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0 .~0.0 .~0.0] ~[.~0.0 .~0.0 .~0.0] ~[.~0.0 .~0.0 .~0.0]]])
  =/  assay-zeros-3x3-6r  (zeros:la meta.input-ones-3x3-6r)
  %+  is-equal
    canon-zeros-3x3-6r
  assay-zeros-3x3-6r
::

++  test-zeros-1x1-7r  ^-  tang
  =/  input-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0]]])
  =/  canon-zeros-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0]]])
  =/  assay-zeros-1x1-7r  (zeros:la meta.input-ones-1x1-7r)
  %+  is-equal
    canon-zeros-1x1-7r
  assay-zeros-1x1-7r
::

++  test-zeros-1x2-7r  ^-  tang
  =/  input-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  canon-zeros-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0 .~~~0.0]]])
  =/  assay-zeros-1x2-7r  (zeros:la meta.input-ones-1x2-7r)
  %+  is-equal
    canon-zeros-1x2-7r
  assay-zeros-1x2-7r
::

++  test-zeros-1x3-7r  ^-  tang
  =/  input-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-zeros-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0 .~~~0.0 .~~~0.0]]])
  =/  assay-zeros-1x3-7r  (zeros:la meta.input-ones-1x3-7r)
  %+  is-equal
    canon-zeros-1x3-7r
  assay-zeros-1x3-7r
::

++  test-zeros-2x1-7r  ^-  tang
  =/  input-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-zeros-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0] ~[.~~~0.0]]])
  =/  assay-zeros-2x1-7r  (zeros:la meta.input-ones-2x1-7r)
  %+  is-equal
    canon-zeros-2x1-7r
  assay-zeros-2x1-7r
::

++  test-zeros-2x2-7r  ^-  tang
  =/  input-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-zeros-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0]]])
  =/  assay-zeros-2x2-7r  (zeros:la meta.input-ones-2x2-7r)
  %+  is-equal
    canon-zeros-2x2-7r
  assay-zeros-2x2-7r
::

++  test-zeros-2x3-7r  ^-  tang
  =/  input-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-zeros-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0 .~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0 .~~~0.0]]])
  =/  assay-zeros-2x3-7r  (zeros:la meta.input-ones-2x3-7r)
  %+  is-equal
    canon-zeros-2x3-7r
  assay-zeros-2x3-7r
::

++  test-zeros-3x1-7r  ^-  tang
  =/  input-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-zeros-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0] ~[.~~~0.0] ~[.~~~0.0]]])
  =/  assay-zeros-3x1-7r  (zeros:la meta.input-ones-3x1-7r)
  %+  is-equal
    canon-zeros-3x1-7r
  assay-zeros-3x1-7r
::

++  test-zeros-3x2-7r  ^-  tang
  =/  input-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-zeros-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0]]])
  =/  assay-zeros-3x2-7r  (zeros:la meta.input-ones-3x2-7r)
  %+  is-equal
    canon-zeros-3x2-7r
  assay-zeros-3x2-7r
::

++  test-zeros-3x3-7r  ^-  tang
  =/  input-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-zeros-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0 .~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0 .~~~0.0]]])
  =/  assay-zeros-3x3-7r  (zeros:la meta.input-ones-3x3-7r)
  %+  is-equal
    canon-zeros-3x3-7r
  assay-zeros-3x3-7r
::

++  test-zeros-1x1-3u  ^-  tang
  =/  input-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-zeros-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-zeros-1x1-3u  (zeros:la meta.input-ones-1x1-3u)
  %+  is-equal
    canon-zeros-1x1-3u
  assay-zeros-1x1-3u
::

++  test-zeros-1x2-3u  ^-  tang
  =/  input-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-zeros-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint prec=~] baum=~[~[0 0]]])
  =/  assay-zeros-1x2-3u  (zeros:la meta.input-ones-1x2-3u)
  %+  is-equal
    canon-zeros-1x2-3u
  assay-zeros-1x2-3u
::

++  test-zeros-1x3-3u  ^-  tang
  =/  input-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  canon-zeros-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint prec=~] baum=~[~[0 0 0]]])
  =/  assay-zeros-1x3-3u  (zeros:la meta.input-ones-1x3-3u)
  %+  is-equal
    canon-zeros-1x3-3u
  assay-zeros-1x3-3u
::

++  test-zeros-2x1-3u  ^-  tang
  =/  input-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  canon-zeros-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  assay-zeros-2x1-3u  (zeros:la meta.input-ones-2x1-3u)
  %+  is-equal
    canon-zeros-2x1-3u
  assay-zeros-2x1-3u
::

++  test-zeros-2x2-3u  ^-  tang
  =/  input-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-zeros-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint prec=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-zeros-2x2-3u  (zeros:la meta.input-ones-2x2-3u)
  %+  is-equal
    canon-zeros-2x2-3u
  assay-zeros-2x2-3u
::

++  test-zeros-2x3-3u  ^-  tang
  =/  input-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-zeros-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint prec=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-zeros-2x3-3u  (zeros:la meta.input-ones-2x3-3u)
  %+  is-equal
    canon-zeros-2x3-3u
  assay-zeros-2x3-3u
::

++  test-zeros-3x1-3u  ^-  tang
  =/  input-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-zeros-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint prec=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-zeros-3x1-3u  (zeros:la meta.input-ones-3x1-3u)
  %+  is-equal
    canon-zeros-3x1-3u
  assay-zeros-3x1-3u
::

++  test-zeros-3x2-3u  ^-  tang
  =/  input-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-zeros-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint prec=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-zeros-3x2-3u  (zeros:la meta.input-ones-3x2-3u)
  %+  is-equal
    canon-zeros-3x2-3u
  assay-zeros-3x2-3u
::

++  test-zeros-3x3-3u  ^-  tang
  =/  input-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-zeros-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint prec=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-zeros-3x3-3u  (zeros:la meta.input-ones-3x3-3u)
  %+  is-equal
    canon-zeros-3x3-3u
  assay-zeros-3x3-3u
::

++  test-zeros-1x1-4u  ^-  tang
  =/  input-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-zeros-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-zeros-1x1-4u  (zeros:la meta.input-ones-1x1-4u)
  %+  is-equal
    canon-zeros-1x1-4u
  assay-zeros-1x1-4u
::

++  test-zeros-1x2-4u  ^-  tang
  =/  input-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-zeros-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint prec=~] baum=~[~[0 0]]])
  =/  assay-zeros-1x2-4u  (zeros:la meta.input-ones-1x2-4u)
  %+  is-equal
    canon-zeros-1x2-4u
  assay-zeros-1x2-4u
::

++  test-zeros-1x3-4u  ^-  tang
  =/  input-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  canon-zeros-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint prec=~] baum=~[~[0 0 0]]])
  =/  assay-zeros-1x3-4u  (zeros:la meta.input-ones-1x3-4u)
  %+  is-equal
    canon-zeros-1x3-4u
  assay-zeros-1x3-4u
::

++  test-zeros-2x1-4u  ^-  tang
  =/  input-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  canon-zeros-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  assay-zeros-2x1-4u  (zeros:la meta.input-ones-2x1-4u)
  %+  is-equal
    canon-zeros-2x1-4u
  assay-zeros-2x1-4u
::

++  test-zeros-2x2-4u  ^-  tang
  =/  input-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-zeros-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint prec=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-zeros-2x2-4u  (zeros:la meta.input-ones-2x2-4u)
  %+  is-equal
    canon-zeros-2x2-4u
  assay-zeros-2x2-4u
::

++  test-zeros-2x3-4u  ^-  tang
  =/  input-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-zeros-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint prec=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-zeros-2x3-4u  (zeros:la meta.input-ones-2x3-4u)
  %+  is-equal
    canon-zeros-2x3-4u
  assay-zeros-2x3-4u
::

++  test-zeros-3x1-4u  ^-  tang
  =/  input-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-zeros-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint prec=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-zeros-3x1-4u  (zeros:la meta.input-ones-3x1-4u)
  %+  is-equal
    canon-zeros-3x1-4u
  assay-zeros-3x1-4u
::

++  test-zeros-3x2-4u  ^-  tang
  =/  input-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-zeros-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint prec=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-zeros-3x2-4u  (zeros:la meta.input-ones-3x2-4u)
  %+  is-equal
    canon-zeros-3x2-4u
  assay-zeros-3x2-4u
::

++  test-zeros-3x3-4u  ^-  tang
  =/  input-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-zeros-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint prec=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-zeros-3x3-4u  (zeros:la meta.input-ones-3x3-4u)
  %+  is-equal
    canon-zeros-3x3-4u
  assay-zeros-3x3-4u
::

++  test-zeros-1x1-5u  ^-  tang
  =/  input-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-zeros-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-zeros-1x1-5u  (zeros:la meta.input-ones-1x1-5u)
  %+  is-equal
    canon-zeros-1x1-5u
  assay-zeros-1x1-5u
::

++  test-zeros-1x2-5u  ^-  tang
  =/  input-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-zeros-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint prec=~] baum=~[~[0 0]]])
  =/  assay-zeros-1x2-5u  (zeros:la meta.input-ones-1x2-5u)
  %+  is-equal
    canon-zeros-1x2-5u
  assay-zeros-1x2-5u
::

++  test-zeros-1x3-5u  ^-  tang
  =/  input-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  canon-zeros-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint prec=~] baum=~[~[0 0 0]]])
  =/  assay-zeros-1x3-5u  (zeros:la meta.input-ones-1x3-5u)
  %+  is-equal
    canon-zeros-1x3-5u
  assay-zeros-1x3-5u
::

++  test-zeros-2x1-5u  ^-  tang
  =/  input-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  canon-zeros-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  assay-zeros-2x1-5u  (zeros:la meta.input-ones-2x1-5u)
  %+  is-equal
    canon-zeros-2x1-5u
  assay-zeros-2x1-5u
::

++  test-zeros-2x2-5u  ^-  tang
  =/  input-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-zeros-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint prec=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-zeros-2x2-5u  (zeros:la meta.input-ones-2x2-5u)
  %+  is-equal
    canon-zeros-2x2-5u
  assay-zeros-2x2-5u
::

++  test-zeros-2x3-5u  ^-  tang
  =/  input-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-zeros-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint prec=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-zeros-2x3-5u  (zeros:la meta.input-ones-2x3-5u)
  %+  is-equal
    canon-zeros-2x3-5u
  assay-zeros-2x3-5u
::

++  test-zeros-3x1-5u  ^-  tang
  =/  input-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-zeros-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint prec=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-zeros-3x1-5u  (zeros:la meta.input-ones-3x1-5u)
  %+  is-equal
    canon-zeros-3x1-5u
  assay-zeros-3x1-5u
::

++  test-zeros-3x2-5u  ^-  tang
  =/  input-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-zeros-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint prec=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-zeros-3x2-5u  (zeros:la meta.input-ones-3x2-5u)
  %+  is-equal
    canon-zeros-3x2-5u
  assay-zeros-3x2-5u
::

++  test-zeros-3x3-5u  ^-  tang
  =/  input-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-zeros-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint prec=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-zeros-3x3-5u  (zeros:la meta.input-ones-3x3-5u)
  %+  is-equal
    canon-zeros-3x3-5u
  assay-zeros-3x3-5u
::

++  test-zeros-1x1-6u  ^-  tang
  =/  input-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-zeros-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-zeros-1x1-6u  (zeros:la meta.input-ones-1x1-6u)
  %+  is-equal
    canon-zeros-1x1-6u
  assay-zeros-1x1-6u
::

++  test-zeros-1x2-6u  ^-  tang
  =/  input-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-zeros-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint prec=~] baum=~[~[0 0]]])
  =/  assay-zeros-1x2-6u  (zeros:la meta.input-ones-1x2-6u)
  %+  is-equal
    canon-zeros-1x2-6u
  assay-zeros-1x2-6u
::

++  test-zeros-1x3-6u  ^-  tang
  =/  input-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  canon-zeros-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint prec=~] baum=~[~[0 0 0]]])
  =/  assay-zeros-1x3-6u  (zeros:la meta.input-ones-1x3-6u)
  %+  is-equal
    canon-zeros-1x3-6u
  assay-zeros-1x3-6u
::

++  test-zeros-2x1-6u  ^-  tang
  =/  input-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  canon-zeros-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  assay-zeros-2x1-6u  (zeros:la meta.input-ones-2x1-6u)
  %+  is-equal
    canon-zeros-2x1-6u
  assay-zeros-2x1-6u
::

++  test-zeros-2x2-6u  ^-  tang
  =/  input-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-zeros-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint prec=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-zeros-2x2-6u  (zeros:la meta.input-ones-2x2-6u)
  %+  is-equal
    canon-zeros-2x2-6u
  assay-zeros-2x2-6u
::

++  test-zeros-2x3-6u  ^-  tang
  =/  input-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-zeros-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint prec=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-zeros-2x3-6u  (zeros:la meta.input-ones-2x3-6u)
  %+  is-equal
    canon-zeros-2x3-6u
  assay-zeros-2x3-6u
::

++  test-zeros-3x1-6u  ^-  tang
  =/  input-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-zeros-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint prec=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-zeros-3x1-6u  (zeros:la meta.input-ones-3x1-6u)
  %+  is-equal
    canon-zeros-3x1-6u
  assay-zeros-3x1-6u
::

++  test-zeros-3x2-6u  ^-  tang
  =/  input-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-zeros-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint prec=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-zeros-3x2-6u  (zeros:la meta.input-ones-3x2-6u)
  %+  is-equal
    canon-zeros-3x2-6u
  assay-zeros-3x2-6u
::

++  test-zeros-3x3-6u  ^-  tang
  =/  input-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-zeros-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint prec=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-zeros-3x3-6u  (zeros:la meta.input-ones-3x3-6u)
  %+  is-equal
    canon-zeros-3x3-6u
  assay-zeros-3x3-6u
::

++  test-ones-1x1-4r  ^-  tang
  =/  input-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0]]])
  =/  canon-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0]]])
  =/  assay-ones-1x1-4r  (ones:la meta.input-ones-1x1-4r)
  %+  is-equal
    canon-ones-1x1-4r
  assay-ones-1x1-4r
::

++  test-ones-1x2-4r  ^-  tang
  =/  input-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  canon-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  assay-ones-1x2-4r  (ones:la meta.input-ones-1x2-4r)
  %+  is-equal
    canon-ones-1x2-4r
  assay-ones-1x2-4r
::

++  test-ones-1x3-4r  ^-  tang
  =/  input-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  assay-ones-1x3-4r  (ones:la meta.input-ones-1x3-4r)
  %+  is-equal
    canon-ones-1x3-4r
  assay-ones-1x3-4r
::

++  test-ones-2x1-4r  ^-  tang
  =/  input-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  canon-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  assay-ones-2x1-4r  (ones:la meta.input-ones-2x1-4r)
  %+  is-equal
    canon-ones-2x1-4r
  assay-ones-2x1-4r
::

++  test-ones-2x2-4r  ^-  tang
  =/  input-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  assay-ones-2x2-4r  (ones:la meta.input-ones-2x2-4r)
  %+  is-equal
    canon-ones-2x2-4r
  assay-ones-2x2-4r
::

++  test-ones-2x3-4r  ^-  tang
  =/  input-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  assay-ones-2x3-4r  (ones:la meta.input-ones-2x3-4r)
  %+  is-equal
    canon-ones-2x3-4r
  assay-ones-2x3-4r
::

++  test-ones-3x1-4r  ^-  tang
  =/  input-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  canon-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  assay-ones-3x1-4r  (ones:la meta.input-ones-3x1-4r)
  %+  is-equal
    canon-ones-3x1-4r
  assay-ones-3x1-4r
::

++  test-ones-3x2-4r  ^-  tang
  =/  input-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  assay-ones-3x2-4r  (ones:la meta.input-ones-3x2-4r)
  %+  is-equal
    canon-ones-3x2-4r
  assay-ones-3x2-4r
::

++  test-ones-3x3-4r  ^-  tang
  =/  input-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  assay-ones-3x3-4r  (ones:la meta.input-ones-3x3-4r)
  %+  is-equal
    canon-ones-3x3-4r
  assay-ones-3x3-4r
::

++  test-ones-1x1-5r  ^-  tang
  =/  input-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0]]])
  =/  canon-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0]]])
  =/  assay-ones-1x1-5r  (ones:la meta.input-ones-1x1-5r)
  %+  is-equal
    canon-ones-1x1-5r
  assay-ones-1x1-5r
::

++  test-ones-1x2-5r  ^-  tang
  =/  input-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0]]])
  =/  canon-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0]]])
  =/  assay-ones-1x2-5r  (ones:la meta.input-ones-1x2-5r)
  %+  is-equal
    canon-ones-1x2-5r
  assay-ones-1x2-5r
::

++  test-ones-1x3-5r  ^-  tang
  =/  input-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  canon-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  assay-ones-1x3-5r  (ones:la meta.input-ones-1x3-5r)
  %+  is-equal
    canon-ones-1x3-5r
  assay-ones-1x3-5r
::

++  test-ones-2x1-5r  ^-  tang
  =/  input-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  canon-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  assay-ones-2x1-5r  (ones:la meta.input-ones-2x1-5r)
  %+  is-equal
    canon-ones-2x1-5r
  assay-ones-2x1-5r
::

++  test-ones-2x2-5r  ^-  tang
  =/  input-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  assay-ones-2x2-5r  (ones:la meta.input-ones-2x2-5r)
  %+  is-equal
    canon-ones-2x2-5r
  assay-ones-2x2-5r
::

++  test-ones-2x3-5r  ^-  tang
  =/  input-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  assay-ones-2x3-5r  (ones:la meta.input-ones-2x3-5r)
  %+  is-equal
    canon-ones-2x3-5r
  assay-ones-2x3-5r
::

++  test-ones-3x1-5r  ^-  tang
  =/  input-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  canon-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  assay-ones-3x1-5r  (ones:la meta.input-ones-3x1-5r)
  %+  is-equal
    canon-ones-3x1-5r
  assay-ones-3x1-5r
::

++  test-ones-3x2-5r  ^-  tang
  =/  input-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  assay-ones-3x2-5r  (ones:la meta.input-ones-3x2-5r)
  %+  is-equal
    canon-ones-3x2-5r
  assay-ones-3x2-5r
::

++  test-ones-3x3-5r  ^-  tang
  =/  input-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  assay-ones-3x3-5r  (ones:la meta.input-ones-3x3-5r)
  %+  is-equal
    canon-ones-3x3-5r
  assay-ones-3x3-5r
::

++  test-ones-1x1-6r  ^-  tang
  =/  input-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0]]])
  =/  canon-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0]]])
  =/  assay-ones-1x1-6r  (ones:la meta.input-ones-1x1-6r)
  %+  is-equal
    canon-ones-1x1-6r
  assay-ones-1x1-6r
::

++  test-ones-1x2-6r  ^-  tang
  =/  input-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0]]])
  =/  canon-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0]]])
  =/  assay-ones-1x2-6r  (ones:la meta.input-ones-1x2-6r)
  %+  is-equal
    canon-ones-1x2-6r
  assay-ones-1x2-6r
::

++  test-ones-1x3-6r  ^-  tang
  =/  input-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  assay-ones-1x3-6r  (ones:la meta.input-ones-1x3-6r)
  %+  is-equal
    canon-ones-1x3-6r
  assay-ones-1x3-6r
::

++  test-ones-2x1-6r  ^-  tang
  =/  input-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  canon-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  assay-ones-2x1-6r  (ones:la meta.input-ones-2x1-6r)
  %+  is-equal
    canon-ones-2x1-6r
  assay-ones-2x1-6r
::

++  test-ones-2x2-6r  ^-  tang
  =/  input-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  assay-ones-2x2-6r  (ones:la meta.input-ones-2x2-6r)
  %+  is-equal
    canon-ones-2x2-6r
  assay-ones-2x2-6r
::

++  test-ones-2x3-6r  ^-  tang
  =/  input-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  assay-ones-2x3-6r  (ones:la meta.input-ones-2x3-6r)
  %+  is-equal
    canon-ones-2x3-6r
  assay-ones-2x3-6r
::

++  test-ones-3x1-6r  ^-  tang
  =/  input-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  canon-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  assay-ones-3x1-6r  (ones:la meta.input-ones-3x1-6r)
  %+  is-equal
    canon-ones-3x1-6r
  assay-ones-3x1-6r
::

++  test-ones-3x2-6r  ^-  tang
  =/  input-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  assay-ones-3x2-6r  (ones:la meta.input-ones-3x2-6r)
  %+  is-equal
    canon-ones-3x2-6r
  assay-ones-3x2-6r
::

++  test-ones-3x3-6r  ^-  tang
  =/  input-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  assay-ones-3x3-6r  (ones:la meta.input-ones-3x3-6r)
  %+  is-equal
    canon-ones-3x3-6r
  assay-ones-3x3-6r
::

++  test-ones-1x1-7r  ^-  tang
  =/  input-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0]]])
  =/  canon-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0]]])
  =/  assay-ones-1x1-7r  (ones:la meta.input-ones-1x1-7r)
  %+  is-equal
    canon-ones-1x1-7r
  assay-ones-1x1-7r
::

++  test-ones-1x2-7r  ^-  tang
  =/  input-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  canon-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  assay-ones-1x2-7r  (ones:la meta.input-ones-1x2-7r)
  %+  is-equal
    canon-ones-1x2-7r
  assay-ones-1x2-7r
::

++  test-ones-1x3-7r  ^-  tang
  =/  input-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  assay-ones-1x3-7r  (ones:la meta.input-ones-1x3-7r)
  %+  is-equal
    canon-ones-1x3-7r
  assay-ones-1x3-7r
::

++  test-ones-2x1-7r  ^-  tang
  =/  input-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  assay-ones-2x1-7r  (ones:la meta.input-ones-2x1-7r)
  %+  is-equal
    canon-ones-2x1-7r
  assay-ones-2x1-7r
::

++  test-ones-2x2-7r  ^-  tang
  =/  input-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  assay-ones-2x2-7r  (ones:la meta.input-ones-2x2-7r)
  %+  is-equal
    canon-ones-2x2-7r
  assay-ones-2x2-7r
::

++  test-ones-2x3-7r  ^-  tang
  =/  input-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  assay-ones-2x3-7r  (ones:la meta.input-ones-2x3-7r)
  %+  is-equal
    canon-ones-2x3-7r
  assay-ones-2x3-7r
::

++  test-ones-3x1-7r  ^-  tang
  =/  input-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  assay-ones-3x1-7r  (ones:la meta.input-ones-3x1-7r)
  %+  is-equal
    canon-ones-3x1-7r
  assay-ones-3x1-7r
::

++  test-ones-3x2-7r  ^-  tang
  =/  input-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  assay-ones-3x2-7r  (ones:la meta.input-ones-3x2-7r)
  %+  is-equal
    canon-ones-3x2-7r
  assay-ones-3x2-7r
::

++  test-ones-3x3-7r  ^-  tang
  =/  input-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  assay-ones-3x3-7r  (ones:la meta.input-ones-3x3-7r)
  %+  is-equal
    canon-ones-3x3-7r
  assay-ones-3x3-7r
::

++  test-ones-1x1-3u  ^-  tang
  =/  input-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[1]]])
  =/  assay-ones-1x1-3u  (ones:la meta.input-ones-1x1-3u)
  %+  is-equal
    canon-ones-1x1-3u
  assay-ones-1x1-3u
::

++  test-ones-1x2-3u  ^-  tang
  =/  input-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  assay-ones-1x2-3u  (ones:la meta.input-ones-1x2-3u)
  %+  is-equal
    canon-ones-1x2-3u
  assay-ones-1x2-3u
::

++  test-ones-1x3-3u  ^-  tang
  =/  input-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  canon-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  assay-ones-1x3-3u  (ones:la meta.input-ones-1x3-3u)
  %+  is-equal
    canon-ones-1x3-3u
  assay-ones-1x3-3u
::

++  test-ones-2x1-3u  ^-  tang
  =/  input-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  canon-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  assay-ones-2x1-3u  (ones:la meta.input-ones-2x1-3u)
  %+  is-equal
    canon-ones-2x1-3u
  assay-ones-2x1-3u
::

++  test-ones-2x2-3u  ^-  tang
  =/  input-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-ones-2x2-3u  (ones:la meta.input-ones-2x2-3u)
  %+  is-equal
    canon-ones-2x2-3u
  assay-ones-2x2-3u
::

++  test-ones-2x3-3u  ^-  tang
  =/  input-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-ones-2x3-3u  (ones:la meta.input-ones-2x3-3u)
  %+  is-equal
    canon-ones-2x3-3u
  assay-ones-2x3-3u
::

++  test-ones-3x1-3u  ^-  tang
  =/  input-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-ones-3x1-3u  (ones:la meta.input-ones-3x1-3u)
  %+  is-equal
    canon-ones-3x1-3u
  assay-ones-3x1-3u
::

++  test-ones-3x2-3u  ^-  tang
  =/  input-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-ones-3x2-3u  (ones:la meta.input-ones-3x2-3u)
  %+  is-equal
    canon-ones-3x2-3u
  assay-ones-3x2-3u
::

++  test-ones-3x3-3u  ^-  tang
  =/  input-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-ones-3x3-3u  (ones:la meta.input-ones-3x3-3u)
  %+  is-equal
    canon-ones-3x3-3u
  assay-ones-3x3-3u
::

++  test-ones-1x1-4u  ^-  tang
  =/  input-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[1]]])
  =/  assay-ones-1x1-4u  (ones:la meta.input-ones-1x1-4u)
  %+  is-equal
    canon-ones-1x1-4u
  assay-ones-1x1-4u
::

++  test-ones-1x2-4u  ^-  tang
  =/  input-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  assay-ones-1x2-4u  (ones:la meta.input-ones-1x2-4u)
  %+  is-equal
    canon-ones-1x2-4u
  assay-ones-1x2-4u
::

++  test-ones-1x3-4u  ^-  tang
  =/  input-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  canon-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  assay-ones-1x3-4u  (ones:la meta.input-ones-1x3-4u)
  %+  is-equal
    canon-ones-1x3-4u
  assay-ones-1x3-4u
::

++  test-ones-2x1-4u  ^-  tang
  =/  input-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  canon-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  assay-ones-2x1-4u  (ones:la meta.input-ones-2x1-4u)
  %+  is-equal
    canon-ones-2x1-4u
  assay-ones-2x1-4u
::

++  test-ones-2x2-4u  ^-  tang
  =/  input-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-ones-2x2-4u  (ones:la meta.input-ones-2x2-4u)
  %+  is-equal
    canon-ones-2x2-4u
  assay-ones-2x2-4u
::

++  test-ones-2x3-4u  ^-  tang
  =/  input-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-ones-2x3-4u  (ones:la meta.input-ones-2x3-4u)
  %+  is-equal
    canon-ones-2x3-4u
  assay-ones-2x3-4u
::

++  test-ones-3x1-4u  ^-  tang
  =/  input-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-ones-3x1-4u  (ones:la meta.input-ones-3x1-4u)
  %+  is-equal
    canon-ones-3x1-4u
  assay-ones-3x1-4u
::

++  test-ones-3x2-4u  ^-  tang
  =/  input-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-ones-3x2-4u  (ones:la meta.input-ones-3x2-4u)
  %+  is-equal
    canon-ones-3x2-4u
  assay-ones-3x2-4u
::

++  test-ones-3x3-4u  ^-  tang
  =/  input-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-ones-3x3-4u  (ones:la meta.input-ones-3x3-4u)
  %+  is-equal
    canon-ones-3x3-4u
  assay-ones-3x3-4u
::

++  test-ones-1x1-5u  ^-  tang
  =/  input-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[1]]])
  =/  assay-ones-1x1-5u  (ones:la meta.input-ones-1x1-5u)
  %+  is-equal
    canon-ones-1x1-5u
  assay-ones-1x1-5u
::

++  test-ones-1x2-5u  ^-  tang
  =/  input-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  assay-ones-1x2-5u  (ones:la meta.input-ones-1x2-5u)
  %+  is-equal
    canon-ones-1x2-5u
  assay-ones-1x2-5u
::

++  test-ones-1x3-5u  ^-  tang
  =/  input-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  canon-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  assay-ones-1x3-5u  (ones:la meta.input-ones-1x3-5u)
  %+  is-equal
    canon-ones-1x3-5u
  assay-ones-1x3-5u
::

++  test-ones-2x1-5u  ^-  tang
  =/  input-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  canon-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  assay-ones-2x1-5u  (ones:la meta.input-ones-2x1-5u)
  %+  is-equal
    canon-ones-2x1-5u
  assay-ones-2x1-5u
::

++  test-ones-2x2-5u  ^-  tang
  =/  input-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-ones-2x2-5u  (ones:la meta.input-ones-2x2-5u)
  %+  is-equal
    canon-ones-2x2-5u
  assay-ones-2x2-5u
::

++  test-ones-2x3-5u  ^-  tang
  =/  input-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-ones-2x3-5u  (ones:la meta.input-ones-2x3-5u)
  %+  is-equal
    canon-ones-2x3-5u
  assay-ones-2x3-5u
::

++  test-ones-3x1-5u  ^-  tang
  =/  input-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-ones-3x1-5u  (ones:la meta.input-ones-3x1-5u)
  %+  is-equal
    canon-ones-3x1-5u
  assay-ones-3x1-5u
::

++  test-ones-3x2-5u  ^-  tang
  =/  input-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-ones-3x2-5u  (ones:la meta.input-ones-3x2-5u)
  %+  is-equal
    canon-ones-3x2-5u
  assay-ones-3x2-5u
::

++  test-ones-3x3-5u  ^-  tang
  =/  input-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-ones-3x3-5u  (ones:la meta.input-ones-3x3-5u)
  %+  is-equal
    canon-ones-3x3-5u
  assay-ones-3x3-5u
::

++  test-ones-1x1-6u  ^-  tang
  =/  input-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[1]]])
  =/  canon-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[1]]])
  =/  assay-ones-1x1-6u  (ones:la meta.input-ones-1x1-6u)
  %+  is-equal
    canon-ones-1x1-6u
  assay-ones-1x1-6u
::

++  test-ones-1x2-6u  ^-  tang
  =/  input-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  assay-ones-1x2-6u  (ones:la meta.input-ones-1x2-6u)
  %+  is-equal
    canon-ones-1x2-6u
  assay-ones-1x2-6u
::

++  test-ones-1x3-6u  ^-  tang
  =/  input-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  canon-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint prec=~] baum=~[~[1 1 1]]])
  =/  assay-ones-1x3-6u  (ones:la meta.input-ones-1x3-6u)
  %+  is-equal
    canon-ones-1x3-6u
  assay-ones-1x3-6u
::

++  test-ones-2x1-6u  ^-  tang
  =/  input-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  canon-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint prec=~] baum=~[~[1] ~[1]]])
  =/  assay-ones-2x1-6u  (ones:la meta.input-ones-2x1-6u)
  %+  is-equal
    canon-ones-2x1-6u
  assay-ones-2x1-6u
::

++  test-ones-2x2-6u  ^-  tang
  =/  input-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-ones-2x2-6u  (ones:la meta.input-ones-2x2-6u)
  %+  is-equal
    canon-ones-2x2-6u
  assay-ones-2x2-6u
::

++  test-ones-2x3-6u  ^-  tang
  =/  input-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-ones-2x3-6u  (ones:la meta.input-ones-2x3-6u)
  %+  is-equal
    canon-ones-2x3-6u
  assay-ones-2x3-6u
::

++  test-ones-3x1-6u  ^-  tang
  =/  input-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint prec=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-ones-3x1-6u  (ones:la meta.input-ones-3x1-6u)
  %+  is-equal
    canon-ones-3x1-6u
  assay-ones-3x1-6u
::

++  test-ones-3x2-6u  ^-  tang
  =/  input-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-ones-3x2-6u  (ones:la meta.input-ones-3x2-6u)
  %+  is-equal
    canon-ones-3x2-6u
  assay-ones-3x2-6u
::

++  test-ones-3x3-6u  ^-  tang
  =/  input-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint prec=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-ones-3x3-6u  (ones:la meta.input-ones-3x3-6u)
  %+  is-equal
    canon-ones-3x3-6u
  assay-ones-3x3-6u
::

++  test-linspace-asc  ^-  tang
  =/  canon-linspace-11  (en-ray:la [meta=[shape=~[11] bloq=5 kind=%i754 prec=~] baum=~[0x0 0x3dcc.cccc 0x3e4c.cccc 0x3e99.9999 0x3ecc.cccc 0x3eff.ffff 0x3f19.9999 0x3f33.3332 0x3f4c.cccc 0x3f66.6665 0x3f80.0000]])
  =/  assay-linspace-11  (linspace:la [~[11] 5 %i754 ~] [.0 .1] 11)
  %+  is-equal
    canon-linspace-11
  assay-linspace-11
::

++  test-linspace-des  ^-  tang
  =/  canon-linspace-11  (en-ray:la [meta=[shape=~[11] bloq=5 kind=%i754 prec=~] baum=~[0x3f80.0000 0x3f66.6666 0x3f4c.cccd 0x3f33.3333 0x3f19.999a 0x3f00.0000 0x3ecc.ccce 0x3e99.999c 0x3e4c.ccd0 0x3dcc.ccd8 0x0]])
  =/  assay-linspace-11  (linspace:la [~[11] 5 %i754 ~] [.1 .0] 11)
  %+  is-equal
    canon-linspace-11
  assay-linspace-11
::

++  test-linspace-1-4r  ^-  tang
  =/  canon-linspace-1-4r  (en-ray:la [meta=[shape=~[1] bloq=4 kind=%i754 prec=~] baum=~[.~~0.0]])
  =/  assay-linspace-1-4r  (linspace:la [~[1] 4 %i754 ~] [.~~0.0 .~~1.0] 1)
  %+  is-equal
    canon-linspace-1-4r
  assay-linspace-1-4r
::

++  test-linspace-1-5r  ^-  tang
  =/  canon-linspace-1-5r  (en-ray:la [meta=[shape=~[1] bloq=5 kind=%i754 prec=~] baum=~[.0.0]])
  =/  assay-linspace-1-5r  (linspace:la [~[1] 5 %i754 ~] [.0.0 .1.0] 1)
  %+  is-equal
    canon-linspace-1-5r
  assay-linspace-1-5r
::

++  test-linspace-1-6r  ^-  tang
  =/  canon-linspace-1-6r  (en-ray:la [meta=[shape=~[1] bloq=6 kind=%i754 prec=~] baum=~[.~0.0]])
  =/  assay-linspace-1-6r  (linspace:la [~[1] 6 %i754 ~] [.~0.0 .~1.0] 1)
  %+  is-equal
    canon-linspace-1-6r
  assay-linspace-1-6r
::

++  test-linspace-1-7r  ^-  tang
  =/  canon-linspace-1-7r  (en-ray:la [meta=[shape=~[1] bloq=7 kind=%i754 prec=~] baum=~[.~~~0.0]])
  =/  assay-linspace-1-7r  (linspace:la [~[1] 7 %i754 ~] [.~~~0.0 .~~~1.0] 1)
  %+  is-equal
    canon-linspace-1-7r
  assay-linspace-1-7r
::

:: TODO test 0x0

++  test-range-asc10-4r  ^-  tang
  =/  meta-1x5-4  [~[10] 4 %i754 ~]
  =/  canon-1x5-4  (en-ray:la [meta-1x5-4 ~[.~~0 .~~0.5 .~~1 .~~1.5 .~~2 .~~2.5 .~~3 .~~3.5 .~~4 .~~4.5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-4)
      !>((range:la meta-1x5-4 [.~~0 .~~5] .~~0.5))
  ==
::

++  test-range-asc10-5r  ^-  tang
  =/  meta-1x5-5  [~[10] 5 %i754 ~]
  =/  canon-1x5-5  (en-ray:la [meta-1x5-5 ~[.0 .0.5 .1 .1.5 .2 .2.5 .3 .3.5 .4 .4.5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-5)
      !>((range:la meta-1x5-5 [.0 .5] .0.5))
  ==
::

++  test-range-asc10-6r  ^-  tang
  =/  meta-1x5-6  [~[10] 6 %i754 ~]
  =/  canon-1x5-6  (en-ray:la [meta-1x5-6 ~[.~0 .~0.5 .~1 .~1.5 .~2 .~2.5 .~3 .~3.5 .~4 .~4.5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-6)
      !>((range:la meta-1x5-6 [.~0 .~5] .~0.5))
  ==
::

++  test-range-asc10-7r  ^-  tang
  =/  meta-1x5-7  [~[10] 7 %i754 ~]
  =/  canon-1x5-7  (en-ray:la [meta-1x5-7 ~[.~~~0 .~~~0.5 .~~~1 .~~~1.5 .~~~2 .~~~2.5 .~~~3 .~~~3.5 .~~~4 .~~~4.5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-7)
      !>((range:la meta-1x5-7 [.~~~0 .~~~5] .~~~0.5))
  ==
::

++  test-range-asc11-4r  ^-  tang
  =/  meta-1x5-4  [~[11] 4 %i754 ~]
  =/  canon-1x5-4  (en-ray:la [meta-1x5-4 ~[.~~0 .~~0.5 .~~1 .~~1.5 .~~2 .~~2.5 .~~3 .~~3.5 .~~4 .~~4.5 .~~5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-4)
      !>((range:la meta-1x5-4 [.~~0 .~~5.1] .~~0.5))
  ==
::

++  test-range-asc11-5r  ^-  tang
  =/  meta-1x5-5  [~[11] 5 %i754 ~]
  =/  canon-1x5-5  (en-ray:la [meta-1x5-5 ~[.0 .0.5 .1 .1.5 .2 .2.5 .3 .3.5 .4 .4.5 .5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-5)
      !>((range:la meta-1x5-5 [.0 .5.1] .0.5))
  ==
::

++  test-range-asc11-6r  ^-  tang
  =/  meta-1x5-6  [~[11] 6 %i754 ~]
  =/  canon-1x5-6  (en-ray:la [meta-1x5-6 ~[.~0 .~0.5 .~1 .~1.5 .~2 .~2.5 .~3 .~3.5 .~4 .~4.5 .~5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-6)
      !>((range:la meta-1x5-6 [.~0 .~5.1] .~0.5))
  ==
::

++  test-range-asc11-7r  ^-  tang
  =/  meta-1x5-7  [~[11] 7 %i754 ~]
  =/  canon-1x5-7  (en-ray:la [meta-1x5-7 ~[.~~~0 .~~~0.5 .~~~1 .~~~1.5 .~~~2 .~~~2.5 .~~~3 .~~~3.5 .~~~4 .~~~4.5 .~~~5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-7)
      !>((range:la meta-1x5-7 [.~~~0 .~~~5.1] .~~~0.5))
  ==
::

++  test-range-des10-4r  ^-  tang
  =/  meta-1x5-4  [~[10] 4 %i754 ~]
  =/  canon-1x5-4  (en-ray:la [meta-1x5-4 ~[.~~5 .~~4.5 .~~4 .~~3.5 .~~3 .~~2.5 .~~2 .~~1.5 .~~1 .~~0.5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-4)
      !>((range:la meta-1x5-4 [.~~5 .~~0] .~~-0.5))
  ==
::

++  test-range-des10-5r  ^-  tang
  =/  meta-1x5-5  [~[10] 5 %i754 ~]
  =/  canon-1x5-5  (en-ray:la [meta-1x5-5 ~[.5 .4.5 .4 .3.5 .3 .2.5 .2 .1.5 .1 .0.5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-5)
      !>((range:la meta-1x5-5 [.5 .0] .-0.5))
  ==
::

++  test-range-des10-6r  ^-  tang
  =/  meta-1x5-6  [~[10] 6 %i754 ~]
  =/  canon-1x5-6  (en-ray:la [meta-1x5-6 ~[.~5 .~4.5 .~4 .~3.5 .~3 .~2.5 .~2 .~1.5 .~1 .~0.5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-6)
      !>((range:la meta-1x5-6 [.~5 .~0] .~-0.5))
  ==
::

++  test-range-des10-7r  ^-  tang
  =/  meta-1x5-7  [~[10] 7 %i754 ~]
  =/  canon-1x5-7  (en-ray:la [meta-1x5-7 ~[.~~~5 .~~~4.5 .~~~4 .~~~3.5 .~~~3 .~~~2.5 .~~~2 .~~~1.5 .~~~1 .~~~0.5]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-7)
      !>((range:la meta-1x5-7 [.~~~5 .~~~0] .~~~-0.5))
  ==
::

++  test-range-des11-4r  ^-  tang
  =/  meta-1x5-4  [~[11] 4 %i754 ~]
  =/  canon-1x5-4  (en-ray:la [meta-1x5-4 ~[.~~5 .~~4.5 .~~4 .~~3.5 .~~3 .~~2.5 .~~2 .~~1.5 .~~1 .~~0.5 .~~0]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-4)
      !>((range:la meta-1x5-4 [.~~5 .~~-0.1] .~~-0.5))
  ==
::

++  test-range-des11-5r  ^-  tang
  =/  meta-1x5-5  [~[11] 5 %i754 ~]
  =/  canon-1x5-5  (en-ray:la [meta-1x5-5 ~[.5 .4.5 .4 .3.5 .3 .2.5 .2 .1.5 .1 .0.5 .0]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-5)
      !>((range:la meta-1x5-5 [.5 .-0.1] .-0.5))
  ==
::

++  test-range-des11-6r  ^-  tang
  =/  meta-1x5-6  [~[11] 6 %i754 ~]
  =/  canon-1x5-6  (en-ray:la [meta-1x5-6 ~[.~5 .~4.5 .~4 .~3.5 .~3 .~2.5 .~2 .~1.5 .~1 .~0.5 .~0]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-6)
      !>((range:la meta-1x5-6 [.~5 .~-0.1] .~-0.5))
  ==
::

++  test-range-des11-7r  ^-  tang
  =/  meta-1x5-7  [~[11] 7 %i754 ~]
  =/  canon-1x5-7  (en-ray:la [meta-1x5-7 ~[.~~~5 .~~~4.5 .~~~4 .~~~3.5 .~~~3 .~~~2.5 .~~~2 .~~~1.5 .~~~1 .~~~0.5 .~~~0]])
  ;:  weld
    %+  expect-eq
      !>(canon-1x5-7)
      !>((range:la meta-1x5-7 [.~~~5 .~~~-0.1] .~~~-0.5))
  ==
--
