/-  *lagoon
/+  *test
/+  *lagoon
::::  /tests/lib/lagoon-compare-reduce -- comparisons and reductions
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
::
++  test-lte-1x1-4r  ^-  tang
  =/  input-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0]]])
  =/  jnput-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0]]])
  =/  canon-lte-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0]]])
  =/  assay-lte-1x1-4r  (lte:la input-ones-1x1-4r jnput-ones-1x1-4r)
  %+  is-equal
    canon-lte-1x1-4r
  assay-lte-1x1-4r
::
::
++  test-lte-1x2-4r  ^-  tang
  =/  input-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  canon-lte-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  assay-lte-1x2-4r  (lte:la input-ones-1x2-4r jnput-ones-1x2-4r)
  %+  is-equal
    canon-lte-1x2-4r
  assay-lte-1x2-4r
::
::
++  test-lte-1x3-4r  ^-  tang
  =/  input-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-lte-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  assay-lte-1x3-4r  (lte:la input-ones-1x3-4r jnput-ones-1x3-4r)
  %+  is-equal
    canon-lte-1x3-4r
  assay-lte-1x3-4r
::
::
++  test-lte-2x1-4r  ^-  tang
  =/  input-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  jnput-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  canon-lte-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  assay-lte-2x1-4r  (lte:la input-ones-2x1-4r jnput-ones-2x1-4r)
  %+  is-equal
    canon-lte-2x1-4r
  assay-lte-2x1-4r
::
::
++  test-lte-2x2-4r  ^-  tang
  =/  input-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-lte-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  assay-lte-2x2-4r  (lte:la input-ones-2x2-4r jnput-ones-2x2-4r)
  %+  is-equal
    canon-lte-2x2-4r
  assay-lte-2x2-4r
::
::
++  test-lte-2x3-4r  ^-  tang
  =/  input-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-lte-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  assay-lte-2x3-4r  (lte:la input-ones-2x3-4r jnput-ones-2x3-4r)
  %+  is-equal
    canon-lte-2x3-4r
  assay-lte-2x3-4r
::
::
++  test-lte-3x1-4r  ^-  tang
  =/  input-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  jnput-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  canon-lte-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  assay-lte-3x1-4r  (lte:la input-ones-3x1-4r jnput-ones-3x1-4r)
  %+  is-equal
    canon-lte-3x1-4r
  assay-lte-3x1-4r
::
::
++  test-lte-3x2-4r  ^-  tang
  =/  input-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-lte-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  assay-lte-3x2-4r  (lte:la input-ones-3x2-4r jnput-ones-3x2-4r)
  %+  is-equal
    canon-lte-3x2-4r
  assay-lte-3x2-4r
::
::
++  test-lte-3x3-4r  ^-  tang
  =/  input-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-lte-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  assay-lte-3x3-4r  (lte:la input-ones-3x3-4r jnput-ones-3x3-4r)
  %+  is-equal
    canon-lte-3x3-4r
  assay-lte-3x3-4r
::
::
++  test-lte-1x1-5r  ^-  tang
  =/  input-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0]]])
  =/  jnput-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0]]])
  =/  canon-lte-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0]]])
  =/  assay-lte-1x1-5r  (lte:la input-ones-1x1-5r jnput-ones-1x1-5r)
  %+  is-equal
    canon-lte-1x1-5r
  assay-lte-1x1-5r
::
::
++  test-lte-1x2-5r  ^-  tang
  =/  input-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0]]])
  =/  jnput-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0]]])
  =/  canon-lte-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0]]])
  =/  assay-lte-1x2-5r  (lte:la input-ones-1x2-5r jnput-ones-1x2-5r)
  %+  is-equal
    canon-lte-1x2-5r
  assay-lte-1x2-5r
::
::
++  test-lte-1x3-5r  ^-  tang
  =/  input-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  canon-lte-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  assay-lte-1x3-5r  (lte:la input-ones-1x3-5r jnput-ones-1x3-5r)
  %+  is-equal
    canon-lte-1x3-5r
  assay-lte-1x3-5r
::
::
++  test-lte-2x1-5r  ^-  tang
  =/  input-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  jnput-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  canon-lte-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  assay-lte-2x1-5r  (lte:la input-ones-2x1-5r jnput-ones-2x1-5r)
  %+  is-equal
    canon-lte-2x1-5r
  assay-lte-2x1-5r
::
::
++  test-lte-2x2-5r  ^-  tang
  =/  input-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  jnput-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-lte-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  assay-lte-2x2-5r  (lte:la input-ones-2x2-5r jnput-ones-2x2-5r)
  %+  is-equal
    canon-lte-2x2-5r
  assay-lte-2x2-5r
::
::
++  test-lte-2x3-5r  ^-  tang
  =/  input-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-lte-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  assay-lte-2x3-5r  (lte:la input-ones-2x3-5r jnput-ones-2x3-5r)
  %+  is-equal
    canon-lte-2x3-5r
  assay-lte-2x3-5r
::
::
++  test-lte-3x1-5r  ^-  tang
  =/  input-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  jnput-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  canon-lte-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  assay-lte-3x1-5r  (lte:la input-ones-3x1-5r jnput-ones-3x1-5r)
  %+  is-equal
    canon-lte-3x1-5r
  assay-lte-3x1-5r
::
::
++  test-lte-3x2-5r  ^-  tang
  =/  input-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  jnput-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-lte-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  assay-lte-3x2-5r  (lte:la input-ones-3x2-5r jnput-ones-3x2-5r)
  %+  is-equal
    canon-lte-3x2-5r
  assay-lte-3x2-5r
::
::
++  test-lte-3x3-5r  ^-  tang
  =/  input-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-lte-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  assay-lte-3x3-5r  (lte:la input-ones-3x3-5r jnput-ones-3x3-5r)
  %+  is-equal
    canon-lte-3x3-5r
  assay-lte-3x3-5r
::
::
++  test-lte-1x1-6r  ^-  tang
  =/  input-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0]]])
  =/  jnput-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0]]])
  =/  canon-lte-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0]]])
  =/  assay-lte-1x1-6r  (lte:la input-ones-1x1-6r jnput-ones-1x1-6r)
  %+  is-equal
    canon-lte-1x1-6r
  assay-lte-1x1-6r
::
::
++  test-lte-1x2-6r  ^-  tang
  =/  input-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0]]])
  =/  jnput-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0]]])
  =/  canon-lte-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0]]])
  =/  assay-lte-1x2-6r  (lte:la input-ones-1x2-6r jnput-ones-1x2-6r)
  %+  is-equal
    canon-lte-1x2-6r
  assay-lte-1x2-6r
::
::
++  test-lte-1x3-6r  ^-  tang
  =/  input-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-lte-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  assay-lte-1x3-6r  (lte:la input-ones-1x3-6r jnput-ones-1x3-6r)
  %+  is-equal
    canon-lte-1x3-6r
  assay-lte-1x3-6r
::
::
++  test-lte-2x1-6r  ^-  tang
  =/  input-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  jnput-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  canon-lte-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  assay-lte-2x1-6r  (lte:la input-ones-2x1-6r jnput-ones-2x1-6r)
  %+  is-equal
    canon-lte-2x1-6r
  assay-lte-2x1-6r
::
::
++  test-lte-2x2-6r  ^-  tang
  =/  input-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  jnput-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-lte-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  assay-lte-2x2-6r  (lte:la input-ones-2x2-6r jnput-ones-2x2-6r)
  %+  is-equal
    canon-lte-2x2-6r
  assay-lte-2x2-6r
::
::
++  test-lte-2x3-6r  ^-  tang
  =/  input-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-lte-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  assay-lte-2x3-6r  (lte:la input-ones-2x3-6r jnput-ones-2x3-6r)
  %+  is-equal
    canon-lte-2x3-6r
  assay-lte-2x3-6r
::
::
++  test-lte-3x1-6r  ^-  tang
  =/  input-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  jnput-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  canon-lte-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  assay-lte-3x1-6r  (lte:la input-ones-3x1-6r jnput-ones-3x1-6r)
  %+  is-equal
    canon-lte-3x1-6r
  assay-lte-3x1-6r
::
::
++  test-lte-3x2-6r  ^-  tang
  =/  input-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  jnput-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-lte-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  assay-lte-3x2-6r  (lte:la input-ones-3x2-6r jnput-ones-3x2-6r)
  %+  is-equal
    canon-lte-3x2-6r
  assay-lte-3x2-6r
::
::
++  test-lte-3x3-6r  ^-  tang
  =/  input-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-lte-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  assay-lte-3x3-6r  (lte:la input-ones-3x3-6r jnput-ones-3x3-6r)
  %+  is-equal
    canon-lte-3x3-6r
  assay-lte-3x3-6r
::
::
++  test-lte-1x1-7r  ^-  tang
  =/  input-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0]]])
  =/  jnput-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0]]])
  =/  canon-lte-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0]]])
  =/  assay-lte-1x1-7r  (lte:la input-ones-1x1-7r jnput-ones-1x1-7r)
  %+  is-equal
    canon-lte-1x1-7r
  assay-lte-1x1-7r
::
::
++  test-lte-1x2-7r  ^-  tang
  =/  input-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  canon-lte-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  assay-lte-1x2-7r  (lte:la input-ones-1x2-7r jnput-ones-1x2-7r)
  %+  is-equal
    canon-lte-1x2-7r
  assay-lte-1x2-7r
::
::
++  test-lte-1x3-7r  ^-  tang
  =/  input-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-lte-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  assay-lte-1x3-7r  (lte:la input-ones-1x3-7r jnput-ones-1x3-7r)
  %+  is-equal
    canon-lte-1x3-7r
  assay-lte-1x3-7r
::
::
++  test-lte-2x1-7r  ^-  tang
  =/  input-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  jnput-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-lte-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  assay-lte-2x1-7r  (lte:la input-ones-2x1-7r jnput-ones-2x1-7r)
  %+  is-equal
    canon-lte-2x1-7r
  assay-lte-2x1-7r
::
::
++  test-lte-2x2-7r  ^-  tang
  =/  input-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-lte-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  assay-lte-2x2-7r  (lte:la input-ones-2x2-7r jnput-ones-2x2-7r)
  %+  is-equal
    canon-lte-2x2-7r
  assay-lte-2x2-7r
::
::
++  test-lte-2x3-7r  ^-  tang
  =/  input-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-lte-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  assay-lte-2x3-7r  (lte:la input-ones-2x3-7r jnput-ones-2x3-7r)
  %+  is-equal
    canon-lte-2x3-7r
  assay-lte-2x3-7r
::
::
++  test-lte-3x1-7r  ^-  tang
  =/  input-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  jnput-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-lte-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  assay-lte-3x1-7r  (lte:la input-ones-3x1-7r jnput-ones-3x1-7r)
  %+  is-equal
    canon-lte-3x1-7r
  assay-lte-3x1-7r
::
::
++  test-lte-3x2-7r  ^-  tang
  =/  input-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-lte-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  assay-lte-3x2-7r  (lte:la input-ones-3x2-7r jnput-ones-3x2-7r)
  %+  is-equal
    canon-lte-3x2-7r
  assay-lte-3x2-7r
::
::
++  test-lte-3x3-7r  ^-  tang
  =/  input-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-lte-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  assay-lte-3x3-7r  (lte:la input-ones-3x3-7r jnput-ones-3x3-7r)
  %+  is-equal
    canon-lte-3x3-7r
  assay-lte-3x3-7r
::
::
++  test-lte-1x1-3u  ^-  tang
  =/  input-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-lte-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[1]]])
  =/  assay-lte-1x1-3u  (lte:la input-ones-1x1-3u jnput-ones-1x1-3u)
  %+  is-equal
    canon-lte-1x1-3u
  assay-lte-1x1-3u
::
::
++  test-lte-1x2-3u  ^-  tang
  =/  input-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-lte-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  assay-lte-1x2-3u  (lte:la input-ones-1x2-3u jnput-ones-1x2-3u)
  %+  is-equal
    canon-lte-1x2-3u
  assay-lte-1x2-3u
::
::
++  test-lte-1x3-3u  ^-  tang
  =/  input-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-lte-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  assay-lte-1x3-3u  (lte:la input-ones-1x3-3u jnput-ones-1x3-3u)
  %+  is-equal
    canon-lte-1x3-3u
  assay-lte-1x3-3u
::
::
++  test-lte-2x1-3u  ^-  tang
  =/  input-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-lte-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  assay-lte-2x1-3u  (lte:la input-ones-2x1-3u jnput-ones-2x1-3u)
  %+  is-equal
    canon-lte-2x1-3u
  assay-lte-2x1-3u
::
::
++  test-lte-2x2-3u  ^-  tang
  =/  input-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-lte-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-lte-2x2-3u  (lte:la input-ones-2x2-3u jnput-ones-2x2-3u)
  %+  is-equal
    canon-lte-2x2-3u
  assay-lte-2x2-3u
::
::
++  test-lte-2x3-3u  ^-  tang
  =/  input-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-lte-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-lte-2x3-3u  (lte:la input-ones-2x3-3u jnput-ones-2x3-3u)
  %+  is-equal
    canon-lte-2x3-3u
  assay-lte-2x3-3u
::
::
++  test-lte-3x1-3u  ^-  tang
  =/  input-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-lte-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-lte-3x1-3u  (lte:la input-ones-3x1-3u jnput-ones-3x1-3u)
  %+  is-equal
    canon-lte-3x1-3u
  assay-lte-3x1-3u
::
::
++  test-lte-3x2-3u  ^-  tang
  =/  input-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-lte-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-lte-3x2-3u  (lte:la input-ones-3x2-3u jnput-ones-3x2-3u)
  %+  is-equal
    canon-lte-3x2-3u
  assay-lte-3x2-3u
::
::
++  test-lte-3x3-3u  ^-  tang
  =/  input-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-lte-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-lte-3x3-3u  (lte:la input-ones-3x3-3u jnput-ones-3x3-3u)
  %+  is-equal
    canon-lte-3x3-3u
  assay-lte-3x3-3u
::
::
++  test-lte-1x1-4u  ^-  tang
  =/  input-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-lte-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[1]]])
  =/  assay-lte-1x1-4u  (lte:la input-ones-1x1-4u jnput-ones-1x1-4u)
  %+  is-equal
    canon-lte-1x1-4u
  assay-lte-1x1-4u
::
::
++  test-lte-1x2-4u  ^-  tang
  =/  input-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-lte-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  assay-lte-1x2-4u  (lte:la input-ones-1x2-4u jnput-ones-1x2-4u)
  %+  is-equal
    canon-lte-1x2-4u
  assay-lte-1x2-4u
::
::
++  test-lte-1x3-4u  ^-  tang
  =/  input-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-lte-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  assay-lte-1x3-4u  (lte:la input-ones-1x3-4u jnput-ones-1x3-4u)
  %+  is-equal
    canon-lte-1x3-4u
  assay-lte-1x3-4u
::
::
++  test-lte-2x1-4u  ^-  tang
  =/  input-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-lte-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  assay-lte-2x1-4u  (lte:la input-ones-2x1-4u jnput-ones-2x1-4u)
  %+  is-equal
    canon-lte-2x1-4u
  assay-lte-2x1-4u
::
::
++  test-lte-2x2-4u  ^-  tang
  =/  input-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-lte-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-lte-2x2-4u  (lte:la input-ones-2x2-4u jnput-ones-2x2-4u)
  %+  is-equal
    canon-lte-2x2-4u
  assay-lte-2x2-4u
::
::
++  test-lte-2x3-4u  ^-  tang
  =/  input-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-lte-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-lte-2x3-4u  (lte:la input-ones-2x3-4u jnput-ones-2x3-4u)
  %+  is-equal
    canon-lte-2x3-4u
  assay-lte-2x3-4u
::
::
++  test-lte-3x1-4u  ^-  tang
  =/  input-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-lte-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-lte-3x1-4u  (lte:la input-ones-3x1-4u jnput-ones-3x1-4u)
  %+  is-equal
    canon-lte-3x1-4u
  assay-lte-3x1-4u
::
::
++  test-lte-3x2-4u  ^-  tang
  =/  input-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-lte-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-lte-3x2-4u  (lte:la input-ones-3x2-4u jnput-ones-3x2-4u)
  %+  is-equal
    canon-lte-3x2-4u
  assay-lte-3x2-4u
::
::
++  test-lte-3x3-4u  ^-  tang
  =/  input-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-lte-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-lte-3x3-4u  (lte:la input-ones-3x3-4u jnput-ones-3x3-4u)
  %+  is-equal
    canon-lte-3x3-4u
  assay-lte-3x3-4u
::
::
++  test-lte-1x1-5u  ^-  tang
  =/  input-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-lte-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[1]]])
  =/  assay-lte-1x1-5u  (lte:la input-ones-1x1-5u jnput-ones-1x1-5u)
  %+  is-equal
    canon-lte-1x1-5u
  assay-lte-1x1-5u
::
::
++  test-lte-1x2-5u  ^-  tang
  =/  input-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-lte-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  assay-lte-1x2-5u  (lte:la input-ones-1x2-5u jnput-ones-1x2-5u)
  %+  is-equal
    canon-lte-1x2-5u
  assay-lte-1x2-5u
::
::
++  test-lte-1x3-5u  ^-  tang
  =/  input-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-lte-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  assay-lte-1x3-5u  (lte:la input-ones-1x3-5u jnput-ones-1x3-5u)
  %+  is-equal
    canon-lte-1x3-5u
  assay-lte-1x3-5u
::
::
++  test-lte-2x1-5u  ^-  tang
  =/  input-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-lte-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  assay-lte-2x1-5u  (lte:la input-ones-2x1-5u jnput-ones-2x1-5u)
  %+  is-equal
    canon-lte-2x1-5u
  assay-lte-2x1-5u
::
::
++  test-lte-2x2-5u  ^-  tang
  =/  input-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-lte-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-lte-2x2-5u  (lte:la input-ones-2x2-5u jnput-ones-2x2-5u)
  %+  is-equal
    canon-lte-2x2-5u
  assay-lte-2x2-5u
::
::
++  test-lte-2x3-5u  ^-  tang
  =/  input-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-lte-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-lte-2x3-5u  (lte:la input-ones-2x3-5u jnput-ones-2x3-5u)
  %+  is-equal
    canon-lte-2x3-5u
  assay-lte-2x3-5u
::
::
++  test-lte-3x1-5u  ^-  tang
  =/  input-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-lte-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-lte-3x1-5u  (lte:la input-ones-3x1-5u jnput-ones-3x1-5u)
  %+  is-equal
    canon-lte-3x1-5u
  assay-lte-3x1-5u
::
::
++  test-lte-3x2-5u  ^-  tang
  =/  input-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-lte-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-lte-3x2-5u  (lte:la input-ones-3x2-5u jnput-ones-3x2-5u)
  %+  is-equal
    canon-lte-3x2-5u
  assay-lte-3x2-5u
::
::
++  test-lte-3x3-5u  ^-  tang
  =/  input-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-lte-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-lte-3x3-5u  (lte:la input-ones-3x3-5u jnput-ones-3x3-5u)
  %+  is-equal
    canon-lte-3x3-5u
  assay-lte-3x3-5u
::
::
++  test-lte-1x1-6u  ^-  tang
  =/  input-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-lte-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[1]]])
  =/  assay-lte-1x1-6u  (lte:la input-ones-1x1-6u jnput-ones-1x1-6u)
  %+  is-equal
    canon-lte-1x1-6u
  assay-lte-1x1-6u
::
::
++  test-lte-1x2-6u  ^-  tang
  =/  input-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-lte-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  assay-lte-1x2-6u  (lte:la input-ones-1x2-6u jnput-ones-1x2-6u)
  %+  is-equal
    canon-lte-1x2-6u
  assay-lte-1x2-6u
::
::
++  test-lte-1x3-6u  ^-  tang
  =/  input-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-lte-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  assay-lte-1x3-6u  (lte:la input-ones-1x3-6u jnput-ones-1x3-6u)
  %+  is-equal
    canon-lte-1x3-6u
  assay-lte-1x3-6u
::
::
++  test-lte-2x1-6u  ^-  tang
  =/  input-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-lte-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  assay-lte-2x1-6u  (lte:la input-ones-2x1-6u jnput-ones-2x1-6u)
  %+  is-equal
    canon-lte-2x1-6u
  assay-lte-2x1-6u
::
::
++  test-lte-2x2-6u  ^-  tang
  =/  input-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-lte-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-lte-2x2-6u  (lte:la input-ones-2x2-6u jnput-ones-2x2-6u)
  %+  is-equal
    canon-lte-2x2-6u
  assay-lte-2x2-6u
::
::
++  test-lte-2x3-6u  ^-  tang
  =/  input-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-lte-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-lte-2x3-6u  (lte:la input-ones-2x3-6u jnput-ones-2x3-6u)
  %+  is-equal
    canon-lte-2x3-6u
  assay-lte-2x3-6u
::
::
++  test-lte-3x1-6u  ^-  tang
  =/  input-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-lte-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-lte-3x1-6u  (lte:la input-ones-3x1-6u jnput-ones-3x1-6u)
  %+  is-equal
    canon-lte-3x1-6u
  assay-lte-3x1-6u
::
::
++  test-lte-3x2-6u  ^-  tang
  =/  input-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-lte-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-lte-3x2-6u  (lte:la input-ones-3x2-6u jnput-ones-3x2-6u)
  %+  is-equal
    canon-lte-3x2-6u
  assay-lte-3x2-6u
::
::
++  test-lte-3x3-6u  ^-  tang
  =/  input-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-lte-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-lte-3x3-6u  (lte:la input-ones-3x3-6u jnput-ones-3x3-6u)
  %+  is-equal
    canon-lte-3x3-6u
  assay-lte-3x3-6u
::
::
++  test-lth-1x1-4r  ^-  tang
  =/  input-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0]]])
  =/  jnput-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0]]])
  =/  canon-lth-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0]]])
  =/  assay-lth-1x1-4r  (lth:la input-ones-1x1-4r jnput-ones-1x1-4r)
  %+  is-equal
    canon-lth-1x1-4r
  assay-lth-1x1-4r
::
::
++  test-lth-1x2-4r  ^-  tang
  =/  input-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  canon-lth-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0]]])
  =/  assay-lth-1x2-4r  (lth:la input-ones-1x2-4r jnput-ones-1x2-4r)
  %+  is-equal
    canon-lth-1x2-4r
  assay-lth-1x2-4r
::
::
++  test-lth-1x3-4r  ^-  tang
  =/  input-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-lth-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0 .~~0.0]]])
  =/  assay-lth-1x3-4r  (lth:la input-ones-1x3-4r jnput-ones-1x3-4r)
  %+  is-equal
    canon-lth-1x3-4r
  assay-lth-1x3-4r
::
::
++  test-lth-2x1-4r  ^-  tang
  =/  input-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  jnput-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  canon-lth-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0] ~[.~~0.0]]])
  =/  assay-lth-2x1-4r  (lth:la input-ones-2x1-4r jnput-ones-2x1-4r)
  %+  is-equal
    canon-lth-2x1-4r
  assay-lth-2x1-4r
::
::
++  test-lth-2x2-4r  ^-  tang
  =/  input-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-lth-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0] ~[.~~0.0 .~~0.0]]])
  =/  assay-lth-2x2-4r  (lth:la input-ones-2x2-4r jnput-ones-2x2-4r)
  %+  is-equal
    canon-lth-2x2-4r
  assay-lth-2x2-4r
::
::
++  test-lth-2x3-4r  ^-  tang
  =/  input-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-lth-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0 .~~0.0] ~[.~~0.0 .~~0.0 .~~0.0]]])
  =/  assay-lth-2x3-4r  (lth:la input-ones-2x3-4r jnput-ones-2x3-4r)
  %+  is-equal
    canon-lth-2x3-4r
  assay-lth-2x3-4r
::
::
++  test-lth-3x1-4r  ^-  tang
  =/  input-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  jnput-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  canon-lth-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0] ~[.~~0.0] ~[.~~0.0]]])
  =/  assay-lth-3x1-4r  (lth:la input-ones-3x1-4r jnput-ones-3x1-4r)
  %+  is-equal
    canon-lth-3x1-4r
  assay-lth-3x1-4r
::
::
++  test-lth-3x2-4r  ^-  tang
  =/  input-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-lth-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0] ~[.~~0.0 .~~0.0] ~[.~~0.0 .~~0.0]]])
  =/  assay-lth-3x2-4r  (lth:la input-ones-3x2-4r jnput-ones-3x2-4r)
  %+  is-equal
    canon-lth-3x2-4r
  assay-lth-3x2-4r
::
::
++  test-lth-3x3-4r  ^-  tang
  =/  input-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-lth-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0 .~~0.0] ~[.~~0.0 .~~0.0 .~~0.0] ~[.~~0.0 .~~0.0 .~~0.0]]])
  =/  assay-lth-3x3-4r  (lth:la input-ones-3x3-4r jnput-ones-3x3-4r)
  %+  is-equal
    canon-lth-3x3-4r
  assay-lth-3x3-4r
::
::
++  test-lth-1x1-5r  ^-  tang
  =/  input-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0]]])
  =/  jnput-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0]]])
  =/  canon-lth-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0]]])
  =/  assay-lth-1x1-5r  (lth:la input-ones-1x1-5r jnput-ones-1x1-5r)
  %+  is-equal
    canon-lth-1x1-5r
  assay-lth-1x1-5r
::
::
++  test-lth-1x2-5r  ^-  tang
  =/  input-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0]]])
  =/  jnput-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0]]])
  =/  canon-lth-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0]]])
  =/  assay-lth-1x2-5r  (lth:la input-ones-1x2-5r jnput-ones-1x2-5r)
  %+  is-equal
    canon-lth-1x2-5r
  assay-lth-1x2-5r
::
::
++  test-lth-1x3-5r  ^-  tang
  =/  input-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  canon-lth-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0 .0.0]]])
  =/  assay-lth-1x3-5r  (lth:la input-ones-1x3-5r jnput-ones-1x3-5r)
  %+  is-equal
    canon-lth-1x3-5r
  assay-lth-1x3-5r
::
::
++  test-lth-2x1-5r  ^-  tang
  =/  input-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  jnput-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  canon-lth-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0] ~[.0.0]]])
  =/  assay-lth-2x1-5r  (lth:la input-ones-2x1-5r jnput-ones-2x1-5r)
  %+  is-equal
    canon-lth-2x1-5r
  assay-lth-2x1-5r
::
::
++  test-lth-2x2-5r  ^-  tang
  =/  input-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  jnput-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-lth-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0] ~[.0.0 .0.0]]])
  =/  assay-lth-2x2-5r  (lth:la input-ones-2x2-5r jnput-ones-2x2-5r)
  %+  is-equal
    canon-lth-2x2-5r
  assay-lth-2x2-5r
::
::
++  test-lth-2x3-5r  ^-  tang
  =/  input-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-lth-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0 .0.0] ~[.0.0 .0.0 .0.0]]])
  =/  assay-lth-2x3-5r  (lth:la input-ones-2x3-5r jnput-ones-2x3-5r)
  %+  is-equal
    canon-lth-2x3-5r
  assay-lth-2x3-5r
::
::
++  test-lth-3x1-5r  ^-  tang
  =/  input-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  jnput-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  canon-lth-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0] ~[.0.0] ~[.0.0]]])
  =/  assay-lth-3x1-5r  (lth:la input-ones-3x1-5r jnput-ones-3x1-5r)
  %+  is-equal
    canon-lth-3x1-5r
  assay-lth-3x1-5r
::
::
++  test-lth-3x2-5r  ^-  tang
  =/  input-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  jnput-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-lth-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0] ~[.0.0 .0.0] ~[.0.0 .0.0]]])
  =/  assay-lth-3x2-5r  (lth:la input-ones-3x2-5r jnput-ones-3x2-5r)
  %+  is-equal
    canon-lth-3x2-5r
  assay-lth-3x2-5r
::
::
++  test-lth-3x3-5r  ^-  tang
  =/  input-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-lth-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0 .0.0] ~[.0.0 .0.0 .0.0] ~[.0.0 .0.0 .0.0]]])
  =/  assay-lth-3x3-5r  (lth:la input-ones-3x3-5r jnput-ones-3x3-5r)
  %+  is-equal
    canon-lth-3x3-5r
  assay-lth-3x3-5r
::
::
++  test-lth-1x1-6r  ^-  tang
  =/  input-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0]]])
  =/  jnput-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0]]])
  =/  canon-lth-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0]]])
  =/  assay-lth-1x1-6r  (lth:la input-ones-1x1-6r jnput-ones-1x1-6r)
  %+  is-equal
    canon-lth-1x1-6r
  assay-lth-1x1-6r
::
::
++  test-lth-1x2-6r  ^-  tang
  =/  input-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0]]])
  =/  jnput-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0]]])
  =/  canon-lth-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0]]])
  =/  assay-lth-1x2-6r  (lth:la input-ones-1x2-6r jnput-ones-1x2-6r)
  %+  is-equal
    canon-lth-1x2-6r
  assay-lth-1x2-6r
::
::
++  test-lth-1x3-6r  ^-  tang
  =/  input-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-lth-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0 .~0.0]]])
  =/  assay-lth-1x3-6r  (lth:la input-ones-1x3-6r jnput-ones-1x3-6r)
  %+  is-equal
    canon-lth-1x3-6r
  assay-lth-1x3-6r
::
::
++  test-lth-2x1-6r  ^-  tang
  =/  input-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  jnput-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  canon-lth-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0] ~[.~0.0]]])
  =/  assay-lth-2x1-6r  (lth:la input-ones-2x1-6r jnput-ones-2x1-6r)
  %+  is-equal
    canon-lth-2x1-6r
  assay-lth-2x1-6r
::
::
++  test-lth-2x2-6r  ^-  tang
  =/  input-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  jnput-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-lth-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0] ~[.~0.0 .~0.0]]])
  =/  assay-lth-2x2-6r  (lth:la input-ones-2x2-6r jnput-ones-2x2-6r)
  %+  is-equal
    canon-lth-2x2-6r
  assay-lth-2x2-6r
::
::
++  test-lth-2x3-6r  ^-  tang
  =/  input-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-lth-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0 .~0.0] ~[.~0.0 .~0.0 .~0.0]]])
  =/  assay-lth-2x3-6r  (lth:la input-ones-2x3-6r jnput-ones-2x3-6r)
  %+  is-equal
    canon-lth-2x3-6r
  assay-lth-2x3-6r
::
::
++  test-lth-3x1-6r  ^-  tang
  =/  input-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  jnput-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  canon-lth-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0] ~[.~0.0] ~[.~0.0]]])
  =/  assay-lth-3x1-6r  (lth:la input-ones-3x1-6r jnput-ones-3x1-6r)
  %+  is-equal
    canon-lth-3x1-6r
  assay-lth-3x1-6r
::
::
++  test-lth-3x2-6r  ^-  tang
  =/  input-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  jnput-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-lth-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0] ~[.~0.0 .~0.0] ~[.~0.0 .~0.0]]])
  =/  assay-lth-3x2-6r  (lth:la input-ones-3x2-6r jnput-ones-3x2-6r)
  %+  is-equal
    canon-lth-3x2-6r
  assay-lth-3x2-6r
::
::
++  test-lth-3x3-6r  ^-  tang
  =/  input-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-lth-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0 .~0.0] ~[.~0.0 .~0.0 .~0.0] ~[.~0.0 .~0.0 .~0.0]]])
  =/  assay-lth-3x3-6r  (lth:la input-ones-3x3-6r jnput-ones-3x3-6r)
  %+  is-equal
    canon-lth-3x3-6r
  assay-lth-3x3-6r
::
::
++  test-lth-1x1-7r  ^-  tang
  =/  input-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0]]])
  =/  jnput-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0]]])
  =/  canon-lth-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0]]])
  =/  assay-lth-1x1-7r  (lth:la input-ones-1x1-7r jnput-ones-1x1-7r)
  %+  is-equal
    canon-lth-1x1-7r
  assay-lth-1x1-7r
::
::
++  test-lth-1x2-7r  ^-  tang
  =/  input-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  canon-lth-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0]]])
  =/  assay-lth-1x2-7r  (lth:la input-ones-1x2-7r jnput-ones-1x2-7r)
  %+  is-equal
    canon-lth-1x2-7r
  assay-lth-1x2-7r
::
::
++  test-lth-1x3-7r  ^-  tang
  =/  input-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-lth-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0 .~~~0.0]]])
  =/  assay-lth-1x3-7r  (lth:la input-ones-1x3-7r jnput-ones-1x3-7r)
  %+  is-equal
    canon-lth-1x3-7r
  assay-lth-1x3-7r
::
::
++  test-lth-2x1-7r  ^-  tang
  =/  input-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  jnput-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-lth-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0] ~[.~~~0.0]]])
  =/  assay-lth-2x1-7r  (lth:la input-ones-2x1-7r jnput-ones-2x1-7r)
  %+  is-equal
    canon-lth-2x1-7r
  assay-lth-2x1-7r
::
::
++  test-lth-2x2-7r  ^-  tang
  =/  input-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-lth-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0]]])
  =/  assay-lth-2x2-7r  (lth:la input-ones-2x2-7r jnput-ones-2x2-7r)
  %+  is-equal
    canon-lth-2x2-7r
  assay-lth-2x2-7r
::
::
++  test-lth-2x3-7r  ^-  tang
  =/  input-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-lth-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0 .~~~0.0]]])
  =/  assay-lth-2x3-7r  (lth:la input-ones-2x3-7r jnput-ones-2x3-7r)
  %+  is-equal
    canon-lth-2x3-7r
  assay-lth-2x3-7r
::
::
++  test-lth-3x1-7r  ^-  tang
  =/  input-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  jnput-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-lth-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0] ~[.~~~0.0] ~[.~~~0.0]]])
  =/  assay-lth-3x1-7r  (lth:la input-ones-3x1-7r jnput-ones-3x1-7r)
  %+  is-equal
    canon-lth-3x1-7r
  assay-lth-3x1-7r
::
::
++  test-lth-3x2-7r  ^-  tang
  =/  input-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-lth-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0]]])
  =/  assay-lth-3x2-7r  (lth:la input-ones-3x2-7r jnput-ones-3x2-7r)
  %+  is-equal
    canon-lth-3x2-7r
  assay-lth-3x2-7r
::
::
++  test-lth-3x3-7r  ^-  tang
  =/  input-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-lth-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0 .~~~0.0]]])
  =/  assay-lth-3x3-7r  (lth:la input-ones-3x3-7r jnput-ones-3x3-7r)
  %+  is-equal
    canon-lth-3x3-7r
  assay-lth-3x3-7r
::
::
++  test-lth-1x1-3u  ^-  tang
  =/  input-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-lth-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[0]]])
  =/  assay-lth-1x1-3u  (lth:la input-ones-1x1-3u jnput-ones-1x1-3u)
  %+  is-equal
    canon-lth-1x1-3u
  assay-lth-1x1-3u
::
::
++  test-lth-1x2-3u  ^-  tang
  =/  input-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-lth-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[0 0]]])
  =/  assay-lth-1x2-3u  (lth:la input-ones-1x2-3u jnput-ones-1x2-3u)
  %+  is-equal
    canon-lth-1x2-3u
  assay-lth-1x2-3u
::
::
++  test-lth-1x3-3u  ^-  tang
  =/  input-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-lth-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[0 0 0]]])
  =/  assay-lth-1x3-3u  (lth:la input-ones-1x3-3u jnput-ones-1x3-3u)
  %+  is-equal
    canon-lth-1x3-3u
  assay-lth-1x3-3u
::
::
++  test-lth-2x1-3u  ^-  tang
  =/  input-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-lth-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[0] ~[0]]])
  =/  assay-lth-2x1-3u  (lth:la input-ones-2x1-3u jnput-ones-2x1-3u)
  %+  is-equal
    canon-lth-2x1-3u
  assay-lth-2x1-3u
::
::
++  test-lth-2x2-3u  ^-  tang
  =/  input-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-lth-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-lth-2x2-3u  (lth:la input-ones-2x2-3u jnput-ones-2x2-3u)
  %+  is-equal
    canon-lth-2x2-3u
  assay-lth-2x2-3u
::
::
++  test-lth-2x3-3u  ^-  tang
  =/  input-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-lth-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-lth-2x3-3u  (lth:la input-ones-2x3-3u jnput-ones-2x3-3u)
  %+  is-equal
    canon-lth-2x3-3u
  assay-lth-2x3-3u
::
::
++  test-lth-3x1-3u  ^-  tang
  =/  input-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-lth-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-lth-3x1-3u  (lth:la input-ones-3x1-3u jnput-ones-3x1-3u)
  %+  is-equal
    canon-lth-3x1-3u
  assay-lth-3x1-3u
::
::
++  test-lth-3x2-3u  ^-  tang
  =/  input-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-lth-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-lth-3x2-3u  (lth:la input-ones-3x2-3u jnput-ones-3x2-3u)
  %+  is-equal
    canon-lth-3x2-3u
  assay-lth-3x2-3u
::
::
++  test-lth-3x3-3u  ^-  tang
  =/  input-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-lth-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-lth-3x3-3u  (lth:la input-ones-3x3-3u jnput-ones-3x3-3u)
  %+  is-equal
    canon-lth-3x3-3u
  assay-lth-3x3-3u
::
::
++  test-lth-1x1-4u  ^-  tang
  =/  input-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-lth-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[0]]])
  =/  assay-lth-1x1-4u  (lth:la input-ones-1x1-4u jnput-ones-1x1-4u)
  %+  is-equal
    canon-lth-1x1-4u
  assay-lth-1x1-4u
::
::
++  test-lth-1x2-4u  ^-  tang
  =/  input-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-lth-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[0 0]]])
  =/  assay-lth-1x2-4u  (lth:la input-ones-1x2-4u jnput-ones-1x2-4u)
  %+  is-equal
    canon-lth-1x2-4u
  assay-lth-1x2-4u
::
::
++  test-lth-1x3-4u  ^-  tang
  =/  input-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-lth-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[0 0 0]]])
  =/  assay-lth-1x3-4u  (lth:la input-ones-1x3-4u jnput-ones-1x3-4u)
  %+  is-equal
    canon-lth-1x3-4u
  assay-lth-1x3-4u
::
::
++  test-lth-2x1-4u  ^-  tang
  =/  input-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-lth-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[0] ~[0]]])
  =/  assay-lth-2x1-4u  (lth:la input-ones-2x1-4u jnput-ones-2x1-4u)
  %+  is-equal
    canon-lth-2x1-4u
  assay-lth-2x1-4u
::
::
++  test-lth-2x2-4u  ^-  tang
  =/  input-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-lth-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-lth-2x2-4u  (lth:la input-ones-2x2-4u jnput-ones-2x2-4u)
  %+  is-equal
    canon-lth-2x2-4u
  assay-lth-2x2-4u
::
::
++  test-lth-2x3-4u  ^-  tang
  =/  input-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-lth-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-lth-2x3-4u  (lth:la input-ones-2x3-4u jnput-ones-2x3-4u)
  %+  is-equal
    canon-lth-2x3-4u
  assay-lth-2x3-4u
::
::
++  test-lth-3x1-4u  ^-  tang
  =/  input-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-lth-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-lth-3x1-4u  (lth:la input-ones-3x1-4u jnput-ones-3x1-4u)
  %+  is-equal
    canon-lth-3x1-4u
  assay-lth-3x1-4u
::
::
++  test-lth-3x2-4u  ^-  tang
  =/  input-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-lth-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-lth-3x2-4u  (lth:la input-ones-3x2-4u jnput-ones-3x2-4u)
  %+  is-equal
    canon-lth-3x2-4u
  assay-lth-3x2-4u
::
::
++  test-lth-3x3-4u  ^-  tang
  =/  input-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-lth-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-lth-3x3-4u  (lth:la input-ones-3x3-4u jnput-ones-3x3-4u)
  %+  is-equal
    canon-lth-3x3-4u
  assay-lth-3x3-4u
::
::
++  test-lth-1x1-5u  ^-  tang
  =/  input-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-lth-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[0]]])
  =/  assay-lth-1x1-5u  (lth:la input-ones-1x1-5u jnput-ones-1x1-5u)
  %+  is-equal
    canon-lth-1x1-5u
  assay-lth-1x1-5u
::
::
++  test-lth-1x2-5u  ^-  tang
  =/  input-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-lth-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[0 0]]])
  =/  assay-lth-1x2-5u  (lth:la input-ones-1x2-5u jnput-ones-1x2-5u)
  %+  is-equal
    canon-lth-1x2-5u
  assay-lth-1x2-5u
::
::
++  test-lth-1x3-5u  ^-  tang
  =/  input-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-lth-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[0 0 0]]])
  =/  assay-lth-1x3-5u  (lth:la input-ones-1x3-5u jnput-ones-1x3-5u)
  %+  is-equal
    canon-lth-1x3-5u
  assay-lth-1x3-5u
::
::
++  test-lth-2x1-5u  ^-  tang
  =/  input-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-lth-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[0] ~[0]]])
  =/  assay-lth-2x1-5u  (lth:la input-ones-2x1-5u jnput-ones-2x1-5u)
  %+  is-equal
    canon-lth-2x1-5u
  assay-lth-2x1-5u
::
::
++  test-lth-2x2-5u  ^-  tang
  =/  input-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-lth-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-lth-2x2-5u  (lth:la input-ones-2x2-5u jnput-ones-2x2-5u)
  %+  is-equal
    canon-lth-2x2-5u
  assay-lth-2x2-5u
::
::
++  test-lth-2x3-5u  ^-  tang
  =/  input-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-lth-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-lth-2x3-5u  (lth:la input-ones-2x3-5u jnput-ones-2x3-5u)
  %+  is-equal
    canon-lth-2x3-5u
  assay-lth-2x3-5u
::
::
++  test-lth-3x1-5u  ^-  tang
  =/  input-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-lth-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-lth-3x1-5u  (lth:la input-ones-3x1-5u jnput-ones-3x1-5u)
  %+  is-equal
    canon-lth-3x1-5u
  assay-lth-3x1-5u
::
::
++  test-lth-3x2-5u  ^-  tang
  =/  input-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-lth-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-lth-3x2-5u  (lth:la input-ones-3x2-5u jnput-ones-3x2-5u)
  %+  is-equal
    canon-lth-3x2-5u
  assay-lth-3x2-5u
::
::
++  test-lth-3x3-5u  ^-  tang
  =/  input-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-lth-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-lth-3x3-5u  (lth:la input-ones-3x3-5u jnput-ones-3x3-5u)
  %+  is-equal
    canon-lth-3x3-5u
  assay-lth-3x3-5u
::
::
++  test-lth-1x1-6u  ^-  tang
  =/  input-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-lth-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[0]]])
  =/  assay-lth-1x1-6u  (lth:la input-ones-1x1-6u jnput-ones-1x1-6u)
  %+  is-equal
    canon-lth-1x1-6u
  assay-lth-1x1-6u
::
::
++  test-lth-1x2-6u  ^-  tang
  =/  input-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-lth-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[0 0]]])
  =/  assay-lth-1x2-6u  (lth:la input-ones-1x2-6u jnput-ones-1x2-6u)
  %+  is-equal
    canon-lth-1x2-6u
  assay-lth-1x2-6u
::
::
++  test-lth-1x3-6u  ^-  tang
  =/  input-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-lth-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[0 0 0]]])
  =/  assay-lth-1x3-6u  (lth:la input-ones-1x3-6u jnput-ones-1x3-6u)
  %+  is-equal
    canon-lth-1x3-6u
  assay-lth-1x3-6u
::
::
++  test-lth-2x1-6u  ^-  tang
  =/  input-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-lth-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[0] ~[0]]])
  =/  assay-lth-2x1-6u  (lth:la input-ones-2x1-6u jnput-ones-2x1-6u)
  %+  is-equal
    canon-lth-2x1-6u
  assay-lth-2x1-6u
::
::
++  test-lth-2x2-6u  ^-  tang
  =/  input-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-lth-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-lth-2x2-6u  (lth:la input-ones-2x2-6u jnput-ones-2x2-6u)
  %+  is-equal
    canon-lth-2x2-6u
  assay-lth-2x2-6u
::
::
++  test-lth-2x3-6u  ^-  tang
  =/  input-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-lth-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-lth-2x3-6u  (lth:la input-ones-2x3-6u jnput-ones-2x3-6u)
  %+  is-equal
    canon-lth-2x3-6u
  assay-lth-2x3-6u
::
::
++  test-lth-3x1-6u  ^-  tang
  =/  input-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-lth-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-lth-3x1-6u  (lth:la input-ones-3x1-6u jnput-ones-3x1-6u)
  %+  is-equal
    canon-lth-3x1-6u
  assay-lth-3x1-6u
::
::
++  test-lth-3x2-6u  ^-  tang
  =/  input-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-lth-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-lth-3x2-6u  (lth:la input-ones-3x2-6u jnput-ones-3x2-6u)
  %+  is-equal
    canon-lth-3x2-6u
  assay-lth-3x2-6u
::
::
++  test-lth-3x3-6u  ^-  tang
  =/  input-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-lth-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-lth-3x3-6u  (lth:la input-ones-3x3-6u jnput-ones-3x3-6u)
  %+  is-equal
    canon-lth-3x3-6u
  assay-lth-3x3-6u
::
::
++  test-gte-1x1-4r  ^-  tang
  =/  input-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0]]])
  =/  jnput-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0]]])
  =/  canon-gte-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0]]])
  =/  assay-gte-1x1-4r  (gte:la input-ones-1x1-4r jnput-ones-1x1-4r)
  %+  is-equal
    canon-gte-1x1-4r
  assay-gte-1x1-4r
::
::
++  test-gte-1x2-4r  ^-  tang
  =/  input-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  canon-gte-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  assay-gte-1x2-4r  (gte:la input-ones-1x2-4r jnput-ones-1x2-4r)
  %+  is-equal
    canon-gte-1x2-4r
  assay-gte-1x2-4r
::
::
++  test-gte-1x3-4r  ^-  tang
  =/  input-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-gte-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  assay-gte-1x3-4r  (gte:la input-ones-1x3-4r jnput-ones-1x3-4r)
  %+  is-equal
    canon-gte-1x3-4r
  assay-gte-1x3-4r
::
::
++  test-gte-2x1-4r  ^-  tang
  =/  input-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  jnput-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  canon-gte-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  assay-gte-2x1-4r  (gte:la input-ones-2x1-4r jnput-ones-2x1-4r)
  %+  is-equal
    canon-gte-2x1-4r
  assay-gte-2x1-4r
::
::
++  test-gte-2x2-4r  ^-  tang
  =/  input-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-gte-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  assay-gte-2x2-4r  (gte:la input-ones-2x2-4r jnput-ones-2x2-4r)
  %+  is-equal
    canon-gte-2x2-4r
  assay-gte-2x2-4r
::
::
++  test-gte-2x3-4r  ^-  tang
  =/  input-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-gte-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  assay-gte-2x3-4r  (gte:la input-ones-2x3-4r jnput-ones-2x3-4r)
  %+  is-equal
    canon-gte-2x3-4r
  assay-gte-2x3-4r
::
::
++  test-gte-3x1-4r  ^-  tang
  =/  input-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  jnput-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  canon-gte-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  assay-gte-3x1-4r  (gte:la input-ones-3x1-4r jnput-ones-3x1-4r)
  %+  is-equal
    canon-gte-3x1-4r
  assay-gte-3x1-4r
::
::
++  test-gte-3x2-4r  ^-  tang
  =/  input-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-gte-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  assay-gte-3x2-4r  (gte:la input-ones-3x2-4r jnput-ones-3x2-4r)
  %+  is-equal
    canon-gte-3x2-4r
  assay-gte-3x2-4r
::
::
++  test-gte-3x3-4r  ^-  tang
  =/  input-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-gte-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  assay-gte-3x3-4r  (gte:la input-ones-3x3-4r jnput-ones-3x3-4r)
  %+  is-equal
    canon-gte-3x3-4r
  assay-gte-3x3-4r
::
::
++  test-gte-1x1-5r  ^-  tang
  =/  input-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0]]])
  =/  jnput-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0]]])
  =/  canon-gte-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0]]])
  =/  assay-gte-1x1-5r  (gte:la input-ones-1x1-5r jnput-ones-1x1-5r)
  %+  is-equal
    canon-gte-1x1-5r
  assay-gte-1x1-5r
::
::
++  test-gte-1x2-5r  ^-  tang
  =/  input-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0]]])
  =/  jnput-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0]]])
  =/  canon-gte-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0]]])
  =/  assay-gte-1x2-5r  (gte:la input-ones-1x2-5r jnput-ones-1x2-5r)
  %+  is-equal
    canon-gte-1x2-5r
  assay-gte-1x2-5r
::
::
++  test-gte-1x3-5r  ^-  tang
  =/  input-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  canon-gte-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  assay-gte-1x3-5r  (gte:la input-ones-1x3-5r jnput-ones-1x3-5r)
  %+  is-equal
    canon-gte-1x3-5r
  assay-gte-1x3-5r
::
::
++  test-gte-2x1-5r  ^-  tang
  =/  input-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  jnput-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  canon-gte-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  assay-gte-2x1-5r  (gte:la input-ones-2x1-5r jnput-ones-2x1-5r)
  %+  is-equal
    canon-gte-2x1-5r
  assay-gte-2x1-5r
::
::
++  test-gte-2x2-5r  ^-  tang
  =/  input-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  jnput-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-gte-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  assay-gte-2x2-5r  (gte:la input-ones-2x2-5r jnput-ones-2x2-5r)
  %+  is-equal
    canon-gte-2x2-5r
  assay-gte-2x2-5r
::
::
++  test-gte-2x3-5r  ^-  tang
  =/  input-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-gte-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  assay-gte-2x3-5r  (gte:la input-ones-2x3-5r jnput-ones-2x3-5r)
  %+  is-equal
    canon-gte-2x3-5r
  assay-gte-2x3-5r
::
::
++  test-gte-3x1-5r  ^-  tang
  =/  input-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  jnput-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  canon-gte-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  assay-gte-3x1-5r  (gte:la input-ones-3x1-5r jnput-ones-3x1-5r)
  %+  is-equal
    canon-gte-3x1-5r
  assay-gte-3x1-5r
::
::
++  test-gte-3x2-5r  ^-  tang
  =/  input-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  jnput-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-gte-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  assay-gte-3x2-5r  (gte:la input-ones-3x2-5r jnput-ones-3x2-5r)
  %+  is-equal
    canon-gte-3x2-5r
  assay-gte-3x2-5r
::
::
++  test-gte-3x3-5r  ^-  tang
  =/  input-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-gte-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  assay-gte-3x3-5r  (gte:la input-ones-3x3-5r jnput-ones-3x3-5r)
  %+  is-equal
    canon-gte-3x3-5r
  assay-gte-3x3-5r
::
::
++  test-gte-1x1-6r  ^-  tang
  =/  input-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0]]])
  =/  jnput-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0]]])
  =/  canon-gte-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0]]])
  =/  assay-gte-1x1-6r  (gte:la input-ones-1x1-6r jnput-ones-1x1-6r)
  %+  is-equal
    canon-gte-1x1-6r
  assay-gte-1x1-6r
::
::
++  test-gte-1x2-6r  ^-  tang
  =/  input-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0]]])
  =/  jnput-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0]]])
  =/  canon-gte-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0]]])
  =/  assay-gte-1x2-6r  (gte:la input-ones-1x2-6r jnput-ones-1x2-6r)
  %+  is-equal
    canon-gte-1x2-6r
  assay-gte-1x2-6r
::
::
++  test-gte-1x3-6r  ^-  tang
  =/  input-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-gte-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  assay-gte-1x3-6r  (gte:la input-ones-1x3-6r jnput-ones-1x3-6r)
  %+  is-equal
    canon-gte-1x3-6r
  assay-gte-1x3-6r
::
::
++  test-gte-2x1-6r  ^-  tang
  =/  input-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  jnput-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  canon-gte-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  assay-gte-2x1-6r  (gte:la input-ones-2x1-6r jnput-ones-2x1-6r)
  %+  is-equal
    canon-gte-2x1-6r
  assay-gte-2x1-6r
::
::
++  test-gte-2x2-6r  ^-  tang
  =/  input-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  jnput-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-gte-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  assay-gte-2x2-6r  (gte:la input-ones-2x2-6r jnput-ones-2x2-6r)
  %+  is-equal
    canon-gte-2x2-6r
  assay-gte-2x2-6r
::
::
++  test-gte-2x3-6r  ^-  tang
  =/  input-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-gte-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  assay-gte-2x3-6r  (gte:la input-ones-2x3-6r jnput-ones-2x3-6r)
  %+  is-equal
    canon-gte-2x3-6r
  assay-gte-2x3-6r
::
::
++  test-gte-3x1-6r  ^-  tang
  =/  input-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  jnput-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  canon-gte-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  assay-gte-3x1-6r  (gte:la input-ones-3x1-6r jnput-ones-3x1-6r)
  %+  is-equal
    canon-gte-3x1-6r
  assay-gte-3x1-6r
::
::
++  test-gte-3x2-6r  ^-  tang
  =/  input-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  jnput-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-gte-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  assay-gte-3x2-6r  (gte:la input-ones-3x2-6r jnput-ones-3x2-6r)
  %+  is-equal
    canon-gte-3x2-6r
  assay-gte-3x2-6r
::
::
++  test-gte-3x3-6r  ^-  tang
  =/  input-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-gte-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  assay-gte-3x3-6r  (gte:la input-ones-3x3-6r jnput-ones-3x3-6r)
  %+  is-equal
    canon-gte-3x3-6r
  assay-gte-3x3-6r
::
::
++  test-gte-1x1-7r  ^-  tang
  =/  input-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0]]])
  =/  jnput-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0]]])
  =/  canon-gte-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0]]])
  =/  assay-gte-1x1-7r  (gte:la input-ones-1x1-7r jnput-ones-1x1-7r)
  %+  is-equal
    canon-gte-1x1-7r
  assay-gte-1x1-7r
::
::
++  test-gte-1x2-7r  ^-  tang
  =/  input-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  canon-gte-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  assay-gte-1x2-7r  (gte:la input-ones-1x2-7r jnput-ones-1x2-7r)
  %+  is-equal
    canon-gte-1x2-7r
  assay-gte-1x2-7r
::
::
++  test-gte-1x3-7r  ^-  tang
  =/  input-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-gte-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  assay-gte-1x3-7r  (gte:la input-ones-1x3-7r jnput-ones-1x3-7r)
  %+  is-equal
    canon-gte-1x3-7r
  assay-gte-1x3-7r
::
::
++  test-gte-2x1-7r  ^-  tang
  =/  input-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  jnput-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-gte-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  assay-gte-2x1-7r  (gte:la input-ones-2x1-7r jnput-ones-2x1-7r)
  %+  is-equal
    canon-gte-2x1-7r
  assay-gte-2x1-7r
::
::
++  test-gte-2x2-7r  ^-  tang
  =/  input-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-gte-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  assay-gte-2x2-7r  (gte:la input-ones-2x2-7r jnput-ones-2x2-7r)
  %+  is-equal
    canon-gte-2x2-7r
  assay-gte-2x2-7r
::
::
++  test-gte-2x3-7r  ^-  tang
  =/  input-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-gte-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  assay-gte-2x3-7r  (gte:la input-ones-2x3-7r jnput-ones-2x3-7r)
  %+  is-equal
    canon-gte-2x3-7r
  assay-gte-2x3-7r
::
::
++  test-gte-3x1-7r  ^-  tang
  =/  input-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  jnput-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-gte-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  assay-gte-3x1-7r  (gte:la input-ones-3x1-7r jnput-ones-3x1-7r)
  %+  is-equal
    canon-gte-3x1-7r
  assay-gte-3x1-7r
::
::
++  test-gte-3x2-7r  ^-  tang
  =/  input-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-gte-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  assay-gte-3x2-7r  (gte:la input-ones-3x2-7r jnput-ones-3x2-7r)
  %+  is-equal
    canon-gte-3x2-7r
  assay-gte-3x2-7r
::
::
++  test-gte-3x3-7r  ^-  tang
  =/  input-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-gte-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  assay-gte-3x3-7r  (gte:la input-ones-3x3-7r jnput-ones-3x3-7r)
  %+  is-equal
    canon-gte-3x3-7r
  assay-gte-3x3-7r
::
::
++  test-gte-1x1-3u  ^-  tang
  =/  input-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-gte-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[1]]])
  =/  assay-gte-1x1-3u  (gte:la input-ones-1x1-3u jnput-ones-1x1-3u)
  %+  is-equal
    canon-gte-1x1-3u
  assay-gte-1x1-3u
::
::
++  test-gte-1x2-3u  ^-  tang
  =/  input-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-gte-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  assay-gte-1x2-3u  (gte:la input-ones-1x2-3u jnput-ones-1x2-3u)
  %+  is-equal
    canon-gte-1x2-3u
  assay-gte-1x2-3u
::
::
++  test-gte-1x3-3u  ^-  tang
  =/  input-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-gte-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  assay-gte-1x3-3u  (gte:la input-ones-1x3-3u jnput-ones-1x3-3u)
  %+  is-equal
    canon-gte-1x3-3u
  assay-gte-1x3-3u
::
::
++  test-gte-2x1-3u  ^-  tang
  =/  input-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-gte-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  assay-gte-2x1-3u  (gte:la input-ones-2x1-3u jnput-ones-2x1-3u)
  %+  is-equal
    canon-gte-2x1-3u
  assay-gte-2x1-3u
::
::
++  test-gte-2x2-3u  ^-  tang
  =/  input-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-gte-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-gte-2x2-3u  (gte:la input-ones-2x2-3u jnput-ones-2x2-3u)
  %+  is-equal
    canon-gte-2x2-3u
  assay-gte-2x2-3u
::
::
++  test-gte-2x3-3u  ^-  tang
  =/  input-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-gte-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-gte-2x3-3u  (gte:la input-ones-2x3-3u jnput-ones-2x3-3u)
  %+  is-equal
    canon-gte-2x3-3u
  assay-gte-2x3-3u
::
::
++  test-gte-3x1-3u  ^-  tang
  =/  input-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-gte-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-gte-3x1-3u  (gte:la input-ones-3x1-3u jnput-ones-3x1-3u)
  %+  is-equal
    canon-gte-3x1-3u
  assay-gte-3x1-3u
::
::
++  test-gte-3x2-3u  ^-  tang
  =/  input-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-gte-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-gte-3x2-3u  (gte:la input-ones-3x2-3u jnput-ones-3x2-3u)
  %+  is-equal
    canon-gte-3x2-3u
  assay-gte-3x2-3u
::
::
++  test-gte-3x3-3u  ^-  tang
  =/  input-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-gte-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-gte-3x3-3u  (gte:la input-ones-3x3-3u jnput-ones-3x3-3u)
  %+  is-equal
    canon-gte-3x3-3u
  assay-gte-3x3-3u
::
::
++  test-gte-1x1-4u  ^-  tang
  =/  input-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-gte-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[1]]])
  =/  assay-gte-1x1-4u  (gte:la input-ones-1x1-4u jnput-ones-1x1-4u)
  %+  is-equal
    canon-gte-1x1-4u
  assay-gte-1x1-4u
::
::
++  test-gte-1x2-4u  ^-  tang
  =/  input-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-gte-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  assay-gte-1x2-4u  (gte:la input-ones-1x2-4u jnput-ones-1x2-4u)
  %+  is-equal
    canon-gte-1x2-4u
  assay-gte-1x2-4u
::
::
++  test-gte-1x3-4u  ^-  tang
  =/  input-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-gte-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  assay-gte-1x3-4u  (gte:la input-ones-1x3-4u jnput-ones-1x3-4u)
  %+  is-equal
    canon-gte-1x3-4u
  assay-gte-1x3-4u
::
::
++  test-gte-2x1-4u  ^-  tang
  =/  input-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-gte-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  assay-gte-2x1-4u  (gte:la input-ones-2x1-4u jnput-ones-2x1-4u)
  %+  is-equal
    canon-gte-2x1-4u
  assay-gte-2x1-4u
::
::
++  test-gte-2x2-4u  ^-  tang
  =/  input-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-gte-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-gte-2x2-4u  (gte:la input-ones-2x2-4u jnput-ones-2x2-4u)
  %+  is-equal
    canon-gte-2x2-4u
  assay-gte-2x2-4u
::
::
++  test-gte-2x3-4u  ^-  tang
  =/  input-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-gte-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-gte-2x3-4u  (gte:la input-ones-2x3-4u jnput-ones-2x3-4u)
  %+  is-equal
    canon-gte-2x3-4u
  assay-gte-2x3-4u
::
::
++  test-gte-3x1-4u  ^-  tang
  =/  input-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-gte-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-gte-3x1-4u  (gte:la input-ones-3x1-4u jnput-ones-3x1-4u)
  %+  is-equal
    canon-gte-3x1-4u
  assay-gte-3x1-4u
::
::
++  test-gte-3x2-4u  ^-  tang
  =/  input-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-gte-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-gte-3x2-4u  (gte:la input-ones-3x2-4u jnput-ones-3x2-4u)
  %+  is-equal
    canon-gte-3x2-4u
  assay-gte-3x2-4u
::
::
++  test-gte-3x3-4u  ^-  tang
  =/  input-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-gte-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-gte-3x3-4u  (gte:la input-ones-3x3-4u jnput-ones-3x3-4u)
  %+  is-equal
    canon-gte-3x3-4u
  assay-gte-3x3-4u
::
::
++  test-gte-1x1-5u  ^-  tang
  =/  input-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-gte-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[1]]])
  =/  assay-gte-1x1-5u  (gte:la input-ones-1x1-5u jnput-ones-1x1-5u)
  %+  is-equal
    canon-gte-1x1-5u
  assay-gte-1x1-5u
::
::
++  test-gte-1x2-5u  ^-  tang
  =/  input-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-gte-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  assay-gte-1x2-5u  (gte:la input-ones-1x2-5u jnput-ones-1x2-5u)
  %+  is-equal
    canon-gte-1x2-5u
  assay-gte-1x2-5u
::
::
++  test-gte-1x3-5u  ^-  tang
  =/  input-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-gte-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  assay-gte-1x3-5u  (gte:la input-ones-1x3-5u jnput-ones-1x3-5u)
  %+  is-equal
    canon-gte-1x3-5u
  assay-gte-1x3-5u
::
::
++  test-gte-2x1-5u  ^-  tang
  =/  input-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-gte-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  assay-gte-2x1-5u  (gte:la input-ones-2x1-5u jnput-ones-2x1-5u)
  %+  is-equal
    canon-gte-2x1-5u
  assay-gte-2x1-5u
::
::
++  test-gte-2x2-5u  ^-  tang
  =/  input-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-gte-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-gte-2x2-5u  (gte:la input-ones-2x2-5u jnput-ones-2x2-5u)
  %+  is-equal
    canon-gte-2x2-5u
  assay-gte-2x2-5u
::
::
++  test-gte-2x3-5u  ^-  tang
  =/  input-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-gte-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-gte-2x3-5u  (gte:la input-ones-2x3-5u jnput-ones-2x3-5u)
  %+  is-equal
    canon-gte-2x3-5u
  assay-gte-2x3-5u
::
::
++  test-gte-3x1-5u  ^-  tang
  =/  input-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-gte-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-gte-3x1-5u  (gte:la input-ones-3x1-5u jnput-ones-3x1-5u)
  %+  is-equal
    canon-gte-3x1-5u
  assay-gte-3x1-5u
::
::
++  test-gte-3x2-5u  ^-  tang
  =/  input-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-gte-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-gte-3x2-5u  (gte:la input-ones-3x2-5u jnput-ones-3x2-5u)
  %+  is-equal
    canon-gte-3x2-5u
  assay-gte-3x2-5u
::
::
++  test-gte-3x3-5u  ^-  tang
  =/  input-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-gte-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-gte-3x3-5u  (gte:la input-ones-3x3-5u jnput-ones-3x3-5u)
  %+  is-equal
    canon-gte-3x3-5u
  assay-gte-3x3-5u
::
::
++  test-gte-1x1-6u  ^-  tang
  =/  input-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-gte-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[1]]])
  =/  assay-gte-1x1-6u  (gte:la input-ones-1x1-6u jnput-ones-1x1-6u)
  %+  is-equal
    canon-gte-1x1-6u
  assay-gte-1x1-6u
::
::
++  test-gte-1x2-6u  ^-  tang
  =/  input-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-gte-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  assay-gte-1x2-6u  (gte:la input-ones-1x2-6u jnput-ones-1x2-6u)
  %+  is-equal
    canon-gte-1x2-6u
  assay-gte-1x2-6u
::
::
++  test-gte-1x3-6u  ^-  tang
  =/  input-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-gte-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  assay-gte-1x3-6u  (gte:la input-ones-1x3-6u jnput-ones-1x3-6u)
  %+  is-equal
    canon-gte-1x3-6u
  assay-gte-1x3-6u
::
::
++  test-gte-2x1-6u  ^-  tang
  =/  input-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-gte-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  assay-gte-2x1-6u  (gte:la input-ones-2x1-6u jnput-ones-2x1-6u)
  %+  is-equal
    canon-gte-2x1-6u
  assay-gte-2x1-6u
::
::
++  test-gte-2x2-6u  ^-  tang
  =/  input-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-gte-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  assay-gte-2x2-6u  (gte:la input-ones-2x2-6u jnput-ones-2x2-6u)
  %+  is-equal
    canon-gte-2x2-6u
  assay-gte-2x2-6u
::
::
++  test-gte-2x3-6u  ^-  tang
  =/  input-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-gte-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  assay-gte-2x3-6u  (gte:la input-ones-2x3-6u jnput-ones-2x3-6u)
  %+  is-equal
    canon-gte-2x3-6u
  assay-gte-2x3-6u
::
::
++  test-gte-3x1-6u  ^-  tang
  =/  input-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-gte-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  assay-gte-3x1-6u  (gte:la input-ones-3x1-6u jnput-ones-3x1-6u)
  %+  is-equal
    canon-gte-3x1-6u
  assay-gte-3x1-6u
::
::
++  test-gte-3x2-6u  ^-  tang
  =/  input-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-gte-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  assay-gte-3x2-6u  (gte:la input-ones-3x2-6u jnput-ones-3x2-6u)
  %+  is-equal
    canon-gte-3x2-6u
  assay-gte-3x2-6u
::
::
++  test-gte-3x3-6u  ^-  tang
  =/  input-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-gte-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  assay-gte-3x3-6u  (gte:la input-ones-3x3-6u jnput-ones-3x3-6u)
  %+  is-equal
    canon-gte-3x3-6u
  assay-gte-3x3-6u
::
::
++  test-gth-1x1-4r  ^-  tang
  =/  input-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0]]])
  =/  jnput-ones-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0]]])
  =/  canon-gth-1x1-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0]]])
  =/  assay-gth-1x1-4r  (gth:la input-ones-1x1-4r jnput-ones-1x1-4r)
  %+  is-equal
    canon-gth-1x1-4r
  assay-gth-1x1-4r
::
::
++  test-gth-1x2-4r  ^-  tang
  =/  input-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0]]])
  =/  canon-gth-1x2-4r  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0]]])
  =/  assay-gth-1x2-4r  (gth:la input-ones-1x2-4r jnput-ones-1x2-4r)
  %+  is-equal
    canon-gth-1x2-4r
  assay-gth-1x2-4r
::
::
++  test-gth-1x3-4r  ^-  tang
  =/  input-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-gth-1x3-4r  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0 .~~0.0]]])
  =/  assay-gth-1x3-4r  (gth:la input-ones-1x3-4r jnput-ones-1x3-4r)
  %+  is-equal
    canon-gth-1x3-4r
  assay-gth-1x3-4r
::
::
++  test-gth-2x1-4r  ^-  tang
  =/  input-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  jnput-ones-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0]]])
  =/  canon-gth-2x1-4r  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0] ~[.~~0.0]]])
  =/  assay-gth-2x1-4r  (gth:la input-ones-2x1-4r jnput-ones-2x1-4r)
  %+  is-equal
    canon-gth-2x1-4r
  assay-gth-2x1-4r
::
::
++  test-gth-2x2-4r  ^-  tang
  =/  input-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-gth-2x2-4r  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0] ~[.~~0.0 .~~0.0]]])
  =/  assay-gth-2x2-4r  (gth:la input-ones-2x2-4r jnput-ones-2x2-4r)
  %+  is-equal
    canon-gth-2x2-4r
  assay-gth-2x2-4r
::
::
++  test-gth-2x3-4r  ^-  tang
  =/  input-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-gth-2x3-4r  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0 .~~0.0] ~[.~~0.0 .~~0.0 .~~0.0]]])
  =/  assay-gth-2x3-4r  (gth:la input-ones-2x3-4r jnput-ones-2x3-4r)
  %+  is-equal
    canon-gth-2x3-4r
  assay-gth-2x3-4r
::
::
++  test-gth-3x1-4r  ^-  tang
  =/  input-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  jnput-ones-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0] ~[.~~1.0] ~[.~~1.0]]])
  =/  canon-gth-3x1-4r  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0] ~[.~~0.0] ~[.~~0.0]]])
  =/  assay-gth-3x1-4r  (gth:la input-ones-3x1-4r jnput-ones-3x1-4r)
  %+  is-equal
    canon-gth-3x1-4r
  assay-gth-3x1-4r
::
::
++  test-gth-3x2-4r  ^-  tang
  =/  input-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  jnput-ones-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0] ~[.~~1.0 .~~1.0]]])
  =/  canon-gth-3x2-4r  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0] ~[.~~0.0 .~~0.0] ~[.~~0.0 .~~0.0]]])
  =/  assay-gth-3x2-4r  (gth:la input-ones-3x2-4r jnput-ones-3x2-4r)
  %+  is-equal
    canon-gth-3x2-4r
  assay-gth-3x2-4r
::
::
++  test-gth-3x3-4r  ^-  tang
  =/  input-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  jnput-ones-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0] ~[.~~1.0 .~~1.0 .~~1.0]]])
  =/  canon-gth-3x3-4r  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%i754 tail=~] baum=~[~[.~~0.0 .~~0.0 .~~0.0] ~[.~~0.0 .~~0.0 .~~0.0] ~[.~~0.0 .~~0.0 .~~0.0]]])
  =/  assay-gth-3x3-4r  (gth:la input-ones-3x3-4r jnput-ones-3x3-4r)
  %+  is-equal
    canon-gth-3x3-4r
  assay-gth-3x3-4r
::
::
++  test-gth-1x1-5r  ^-  tang
  =/  input-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0]]])
  =/  jnput-ones-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0]]])
  =/  canon-gth-1x1-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0]]])
  =/  assay-gth-1x1-5r  (gth:la input-ones-1x1-5r jnput-ones-1x1-5r)
  %+  is-equal
    canon-gth-1x1-5r
  assay-gth-1x1-5r
::
::
++  test-gth-1x2-5r  ^-  tang
  =/  input-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0]]])
  =/  jnput-ones-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0]]])
  =/  canon-gth-1x2-5r  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0]]])
  =/  assay-gth-1x2-5r  (gth:la input-ones-1x2-5r jnput-ones-1x2-5r)
  %+  is-equal
    canon-gth-1x2-5r
  assay-gth-1x2-5r
::
::
++  test-gth-1x3-5r  ^-  tang
  =/  input-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  canon-gth-1x3-5r  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0 .0.0]]])
  =/  assay-gth-1x3-5r  (gth:la input-ones-1x3-5r jnput-ones-1x3-5r)
  %+  is-equal
    canon-gth-1x3-5r
  assay-gth-1x3-5r
::
::
++  test-gth-2x1-5r  ^-  tang
  =/  input-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  jnput-ones-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0]]])
  =/  canon-gth-2x1-5r  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0] ~[.0.0]]])
  =/  assay-gth-2x1-5r  (gth:la input-ones-2x1-5r jnput-ones-2x1-5r)
  %+  is-equal
    canon-gth-2x1-5r
  assay-gth-2x1-5r
::
::
++  test-gth-2x2-5r  ^-  tang
  =/  input-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  jnput-ones-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-gth-2x2-5r  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0] ~[.0.0 .0.0]]])
  =/  assay-gth-2x2-5r  (gth:la input-ones-2x2-5r jnput-ones-2x2-5r)
  %+  is-equal
    canon-gth-2x2-5r
  assay-gth-2x2-5r
::
::
++  test-gth-2x3-5r  ^-  tang
  =/  input-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-gth-2x3-5r  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0 .0.0] ~[.0.0 .0.0 .0.0]]])
  =/  assay-gth-2x3-5r  (gth:la input-ones-2x3-5r jnput-ones-2x3-5r)
  %+  is-equal
    canon-gth-2x3-5r
  assay-gth-2x3-5r
::
::
++  test-gth-3x1-5r  ^-  tang
  =/  input-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  jnput-ones-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0] ~[.1.0] ~[.1.0]]])
  =/  canon-gth-3x1-5r  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0] ~[.0.0] ~[.0.0]]])
  =/  assay-gth-3x1-5r  (gth:la input-ones-3x1-5r jnput-ones-3x1-5r)
  %+  is-equal
    canon-gth-3x1-5r
  assay-gth-3x1-5r
::
::
++  test-gth-3x2-5r  ^-  tang
  =/  input-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  jnput-ones-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0] ~[.1.0 .1.0] ~[.1.0 .1.0]]])
  =/  canon-gth-3x2-5r  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0] ~[.0.0 .0.0] ~[.0.0 .0.0]]])
  =/  assay-gth-3x2-5r  (gth:la input-ones-3x2-5r jnput-ones-3x2-5r)
  %+  is-equal
    canon-gth-3x2-5r
  assay-gth-3x2-5r
::
::
++  test-gth-3x3-5r  ^-  tang
  =/  input-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  jnput-ones-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0] ~[.1.0 .1.0 .1.0]]])
  =/  canon-gth-3x3-5r  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0 .0.0] ~[.0.0 .0.0 .0.0] ~[.0.0 .0.0 .0.0]]])
  =/  assay-gth-3x3-5r  (gth:la input-ones-3x3-5r jnput-ones-3x3-5r)
  %+  is-equal
    canon-gth-3x3-5r
  assay-gth-3x3-5r
::
::
++  test-gth-1x1-6r  ^-  tang
  =/  input-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0]]])
  =/  jnput-ones-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0]]])
  =/  canon-gth-1x1-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0]]])
  =/  assay-gth-1x1-6r  (gth:la input-ones-1x1-6r jnput-ones-1x1-6r)
  %+  is-equal
    canon-gth-1x1-6r
  assay-gth-1x1-6r
::
::
++  test-gth-1x2-6r  ^-  tang
  =/  input-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0]]])
  =/  jnput-ones-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0]]])
  =/  canon-gth-1x2-6r  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0]]])
  =/  assay-gth-1x2-6r  (gth:la input-ones-1x2-6r jnput-ones-1x2-6r)
  %+  is-equal
    canon-gth-1x2-6r
  assay-gth-1x2-6r
::
::
++  test-gth-1x3-6r  ^-  tang
  =/  input-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-gth-1x3-6r  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0 .~0.0]]])
  =/  assay-gth-1x3-6r  (gth:la input-ones-1x3-6r jnput-ones-1x3-6r)
  %+  is-equal
    canon-gth-1x3-6r
  assay-gth-1x3-6r
::
::
++  test-gth-2x1-6r  ^-  tang
  =/  input-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  jnput-ones-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0]]])
  =/  canon-gth-2x1-6r  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0] ~[.~0.0]]])
  =/  assay-gth-2x1-6r  (gth:la input-ones-2x1-6r jnput-ones-2x1-6r)
  %+  is-equal
    canon-gth-2x1-6r
  assay-gth-2x1-6r
::
::
++  test-gth-2x2-6r  ^-  tang
  =/  input-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  jnput-ones-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-gth-2x2-6r  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0] ~[.~0.0 .~0.0]]])
  =/  assay-gth-2x2-6r  (gth:la input-ones-2x2-6r jnput-ones-2x2-6r)
  %+  is-equal
    canon-gth-2x2-6r
  assay-gth-2x2-6r
::
::
++  test-gth-2x3-6r  ^-  tang
  =/  input-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-gth-2x3-6r  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0 .~0.0] ~[.~0.0 .~0.0 .~0.0]]])
  =/  assay-gth-2x3-6r  (gth:la input-ones-2x3-6r jnput-ones-2x3-6r)
  %+  is-equal
    canon-gth-2x3-6r
  assay-gth-2x3-6r
::
::
++  test-gth-3x1-6r  ^-  tang
  =/  input-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  jnput-ones-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0] ~[.~1.0] ~[.~1.0]]])
  =/  canon-gth-3x1-6r  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0] ~[.~0.0] ~[.~0.0]]])
  =/  assay-gth-3x1-6r  (gth:la input-ones-3x1-6r jnput-ones-3x1-6r)
  %+  is-equal
    canon-gth-3x1-6r
  assay-gth-3x1-6r
::
::
++  test-gth-3x2-6r  ^-  tang
  =/  input-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  jnput-ones-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0] ~[.~1.0 .~1.0] ~[.~1.0 .~1.0]]])
  =/  canon-gth-3x2-6r  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0] ~[.~0.0 .~0.0] ~[.~0.0 .~0.0]]])
  =/  assay-gth-3x2-6r  (gth:la input-ones-3x2-6r jnput-ones-3x2-6r)
  %+  is-equal
    canon-gth-3x2-6r
  assay-gth-3x2-6r
::
::
++  test-gth-3x3-6r  ^-  tang
  =/  input-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  jnput-ones-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0] ~[.~1.0 .~1.0 .~1.0]]])
  =/  canon-gth-3x3-6r  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%i754 tail=~] baum=~[~[.~0.0 .~0.0 .~0.0] ~[.~0.0 .~0.0 .~0.0] ~[.~0.0 .~0.0 .~0.0]]])
  =/  assay-gth-3x3-6r  (gth:la input-ones-3x3-6r jnput-ones-3x3-6r)
  %+  is-equal
    canon-gth-3x3-6r
  assay-gth-3x3-6r
::
::
++  test-gth-1x1-7r  ^-  tang
  =/  input-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0]]])
  =/  jnput-ones-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0]]])
  =/  canon-gth-1x1-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0]]])
  =/  assay-gth-1x1-7r  (gth:la input-ones-1x1-7r jnput-ones-1x1-7r)
  %+  is-equal
    canon-gth-1x1-7r
  assay-gth-1x1-7r
::
::
++  test-gth-1x2-7r  ^-  tang
  =/  input-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0]]])
  =/  canon-gth-1x2-7r  (en-ray:la [meta=[shape=~[1 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0]]])
  =/  assay-gth-1x2-7r  (gth:la input-ones-1x2-7r jnput-ones-1x2-7r)
  %+  is-equal
    canon-gth-1x2-7r
  assay-gth-1x2-7r
::
::
++  test-gth-1x3-7r  ^-  tang
  =/  input-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-gth-1x3-7r  (en-ray:la [meta=[shape=~[1 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0 .~~~0.0]]])
  =/  assay-gth-1x3-7r  (gth:la input-ones-1x3-7r jnput-ones-1x3-7r)
  %+  is-equal
    canon-gth-1x3-7r
  assay-gth-1x3-7r
::
::
++  test-gth-2x1-7r  ^-  tang
  =/  input-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  jnput-ones-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-gth-2x1-7r  (en-ray:la [meta=[shape=~[2 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0] ~[.~~~0.0]]])
  =/  assay-gth-2x1-7r  (gth:la input-ones-2x1-7r jnput-ones-2x1-7r)
  %+  is-equal
    canon-gth-2x1-7r
  assay-gth-2x1-7r
::
::
++  test-gth-2x2-7r  ^-  tang
  =/  input-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-gth-2x2-7r  (en-ray:la [meta=[shape=~[2 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0]]])
  =/  assay-gth-2x2-7r  (gth:la input-ones-2x2-7r jnput-ones-2x2-7r)
  %+  is-equal
    canon-gth-2x2-7r
  assay-gth-2x2-7r
::
::
++  test-gth-2x3-7r  ^-  tang
  =/  input-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-gth-2x3-7r  (en-ray:la [meta=[shape=~[2 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0 .~~~0.0]]])
  =/  assay-gth-2x3-7r  (gth:la input-ones-2x3-7r jnput-ones-2x3-7r)
  %+  is-equal
    canon-gth-2x3-7r
  assay-gth-2x3-7r
::
::
++  test-gth-3x1-7r  ^-  tang
  =/  input-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  jnput-ones-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0] ~[.~~~1.0] ~[.~~~1.0]]])
  =/  canon-gth-3x1-7r  (en-ray:la [meta=[shape=~[3 1] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0] ~[.~~~0.0] ~[.~~~0.0]]])
  =/  assay-gth-3x1-7r  (gth:la input-ones-3x1-7r jnput-ones-3x1-7r)
  %+  is-equal
    canon-gth-3x1-7r
  assay-gth-3x1-7r
::
::
++  test-gth-3x2-7r  ^-  tang
  =/  input-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  jnput-ones-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0]]])
  =/  canon-gth-3x2-7r  (en-ray:la [meta=[shape=~[3 2] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0]]])
  =/  assay-gth-3x2-7r  (gth:la input-ones-3x2-7r jnput-ones-3x2-7r)
  %+  is-equal
    canon-gth-3x2-7r
  assay-gth-3x2-7r
::
::
++  test-gth-3x3-7r  ^-  tang
  =/  input-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  jnput-ones-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0] ~[.~~~1.0 .~~~1.0 .~~~1.0]]])
  =/  canon-gth-3x3-7r  (en-ray:la [meta=[shape=~[3 3] bloq=7 kind=%i754 tail=~] baum=~[~[.~~~0.0 .~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0 .~~~0.0] ~[.~~~0.0 .~~~0.0 .~~~0.0]]])
  =/  assay-gth-3x3-7r  (gth:la input-ones-3x3-7r jnput-ones-3x3-7r)
  %+  is-equal
    canon-gth-3x3-7r
  assay-gth-3x3-7r
::
::
++  test-gth-1x1-3u  ^-  tang
  =/  input-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-gth-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint tail=~] baum=~[~[0]]])
  =/  assay-gth-1x1-3u  (gth:la input-ones-1x1-3u jnput-ones-1x1-3u)
  %+  is-equal
    canon-gth-1x1-3u
  assay-gth-1x1-3u
::
::
++  test-gth-1x2-3u  ^-  tang
  =/  input-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-gth-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint tail=~] baum=~[~[0 0]]])
  =/  assay-gth-1x2-3u  (gth:la input-ones-1x2-3u jnput-ones-1x2-3u)
  %+  is-equal
    canon-gth-1x2-3u
  assay-gth-1x2-3u
::
::
++  test-gth-1x3-3u  ^-  tang
  =/  input-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-gth-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint tail=~] baum=~[~[0 0 0]]])
  =/  assay-gth-1x3-3u  (gth:la input-ones-1x3-3u jnput-ones-1x3-3u)
  %+  is-equal
    canon-gth-1x3-3u
  assay-gth-1x3-3u
::
::
++  test-gth-2x1-3u  ^-  tang
  =/  input-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-gth-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint tail=~] baum=~[~[0] ~[0]]])
  =/  assay-gth-2x1-3u  (gth:la input-ones-2x1-3u jnput-ones-2x1-3u)
  %+  is-equal
    canon-gth-2x1-3u
  assay-gth-2x1-3u
::
::
++  test-gth-2x2-3u  ^-  tang
  =/  input-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-gth-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint tail=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-gth-2x2-3u  (gth:la input-ones-2x2-3u jnput-ones-2x2-3u)
  %+  is-equal
    canon-gth-2x2-3u
  assay-gth-2x2-3u
::
::
++  test-gth-2x3-3u  ^-  tang
  =/  input-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-gth-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-gth-2x3-3u  (gth:la input-ones-2x3-3u jnput-ones-2x3-3u)
  %+  is-equal
    canon-gth-2x3-3u
  assay-gth-2x3-3u
::
::
++  test-gth-3x1-3u  ^-  tang
  =/  input-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-gth-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint tail=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-gth-3x1-3u  (gth:la input-ones-3x1-3u jnput-ones-3x1-3u)
  %+  is-equal
    canon-gth-3x1-3u
  assay-gth-3x1-3u
::
::
++  test-gth-3x2-3u  ^-  tang
  =/  input-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-gth-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint tail=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-gth-3x2-3u  (gth:la input-ones-3x2-3u jnput-ones-3x2-3u)
  %+  is-equal
    canon-gth-3x2-3u
  assay-gth-3x2-3u
::
::
++  test-gth-3x3-3u  ^-  tang
  =/  input-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-gth-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-gth-3x3-3u  (gth:la input-ones-3x3-3u jnput-ones-3x3-3u)
  %+  is-equal
    canon-gth-3x3-3u
  assay-gth-3x3-3u
::
::
++  test-gth-1x1-4u  ^-  tang
  =/  input-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-gth-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint tail=~] baum=~[~[0]]])
  =/  assay-gth-1x1-4u  (gth:la input-ones-1x1-4u jnput-ones-1x1-4u)
  %+  is-equal
    canon-gth-1x1-4u
  assay-gth-1x1-4u
::
::
++  test-gth-1x2-4u  ^-  tang
  =/  input-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-gth-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint tail=~] baum=~[~[0 0]]])
  =/  assay-gth-1x2-4u  (gth:la input-ones-1x2-4u jnput-ones-1x2-4u)
  %+  is-equal
    canon-gth-1x2-4u
  assay-gth-1x2-4u
::
::
++  test-gth-1x3-4u  ^-  tang
  =/  input-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-gth-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint tail=~] baum=~[~[0 0 0]]])
  =/  assay-gth-1x3-4u  (gth:la input-ones-1x3-4u jnput-ones-1x3-4u)
  %+  is-equal
    canon-gth-1x3-4u
  assay-gth-1x3-4u
::
::
++  test-gth-2x1-4u  ^-  tang
  =/  input-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-gth-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint tail=~] baum=~[~[0] ~[0]]])
  =/  assay-gth-2x1-4u  (gth:la input-ones-2x1-4u jnput-ones-2x1-4u)
  %+  is-equal
    canon-gth-2x1-4u
  assay-gth-2x1-4u
::
::
++  test-gth-2x2-4u  ^-  tang
  =/  input-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-gth-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint tail=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-gth-2x2-4u  (gth:la input-ones-2x2-4u jnput-ones-2x2-4u)
  %+  is-equal
    canon-gth-2x2-4u
  assay-gth-2x2-4u
::
::
++  test-gth-2x3-4u  ^-  tang
  =/  input-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-gth-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-gth-2x3-4u  (gth:la input-ones-2x3-4u jnput-ones-2x3-4u)
  %+  is-equal
    canon-gth-2x3-4u
  assay-gth-2x3-4u
::
::
++  test-gth-3x1-4u  ^-  tang
  =/  input-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-gth-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint tail=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-gth-3x1-4u  (gth:la input-ones-3x1-4u jnput-ones-3x1-4u)
  %+  is-equal
    canon-gth-3x1-4u
  assay-gth-3x1-4u
::
::
++  test-gth-3x2-4u  ^-  tang
  =/  input-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-gth-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint tail=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-gth-3x2-4u  (gth:la input-ones-3x2-4u jnput-ones-3x2-4u)
  %+  is-equal
    canon-gth-3x2-4u
  assay-gth-3x2-4u
::
::
++  test-gth-3x3-4u  ^-  tang
  =/  input-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-gth-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-gth-3x3-4u  (gth:la input-ones-3x3-4u jnput-ones-3x3-4u)
  %+  is-equal
    canon-gth-3x3-4u
  assay-gth-3x3-4u
::
::
++  test-gth-1x1-5u  ^-  tang
  =/  input-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-gth-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint tail=~] baum=~[~[0]]])
  =/  assay-gth-1x1-5u  (gth:la input-ones-1x1-5u jnput-ones-1x1-5u)
  %+  is-equal
    canon-gth-1x1-5u
  assay-gth-1x1-5u
::
::
++  test-gth-1x2-5u  ^-  tang
  =/  input-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-gth-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint tail=~] baum=~[~[0 0]]])
  =/  assay-gth-1x2-5u  (gth:la input-ones-1x2-5u jnput-ones-1x2-5u)
  %+  is-equal
    canon-gth-1x2-5u
  assay-gth-1x2-5u
::
::
++  test-gth-1x3-5u  ^-  tang
  =/  input-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-gth-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint tail=~] baum=~[~[0 0 0]]])
  =/  assay-gth-1x3-5u  (gth:la input-ones-1x3-5u jnput-ones-1x3-5u)
  %+  is-equal
    canon-gth-1x3-5u
  assay-gth-1x3-5u
::
::
++  test-gth-2x1-5u  ^-  tang
  =/  input-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-gth-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint tail=~] baum=~[~[0] ~[0]]])
  =/  assay-gth-2x1-5u  (gth:la input-ones-2x1-5u jnput-ones-2x1-5u)
  %+  is-equal
    canon-gth-2x1-5u
  assay-gth-2x1-5u
::
::
++  test-gth-2x2-5u  ^-  tang
  =/  input-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-gth-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint tail=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-gth-2x2-5u  (gth:la input-ones-2x2-5u jnput-ones-2x2-5u)
  %+  is-equal
    canon-gth-2x2-5u
  assay-gth-2x2-5u
::
::
++  test-gth-2x3-5u  ^-  tang
  =/  input-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-gth-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-gth-2x3-5u  (gth:la input-ones-2x3-5u jnput-ones-2x3-5u)
  %+  is-equal
    canon-gth-2x3-5u
  assay-gth-2x3-5u
::
::
++  test-gth-3x1-5u  ^-  tang
  =/  input-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-gth-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint tail=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-gth-3x1-5u  (gth:la input-ones-3x1-5u jnput-ones-3x1-5u)
  %+  is-equal
    canon-gth-3x1-5u
  assay-gth-3x1-5u
::
::
++  test-gth-3x2-5u  ^-  tang
  =/  input-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-gth-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint tail=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-gth-3x2-5u  (gth:la input-ones-3x2-5u jnput-ones-3x2-5u)
  %+  is-equal
    canon-gth-3x2-5u
  assay-gth-3x2-5u
::
::
++  test-gth-3x3-5u  ^-  tang
  =/  input-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-gth-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-gth-3x3-5u  (gth:la input-ones-3x3-5u jnput-ones-3x3-5u)
  %+  is-equal
    canon-gth-3x3-5u
  assay-gth-3x3-5u
::
::
++  test-gth-1x1-6u  ^-  tang
  =/  input-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[1]]])
  =/  jnput-ones-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[1]]])
  =/  canon-gth-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint tail=~] baum=~[~[0]]])
  =/  assay-gth-1x1-6u  (gth:la input-ones-1x1-6u jnput-ones-1x1-6u)
  %+  is-equal
    canon-gth-1x1-6u
  assay-gth-1x1-6u
::
::
++  test-gth-1x2-6u  ^-  tang
  =/  input-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  jnput-ones-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1]]])
  =/  canon-gth-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint tail=~] baum=~[~[0 0]]])
  =/  assay-gth-1x2-6u  (gth:la input-ones-1x2-6u jnput-ones-1x2-6u)
  %+  is-equal
    canon-gth-1x2-6u
  assay-gth-1x2-6u
::
::
++  test-gth-1x3-6u  ^-  tang
  =/  input-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  jnput-ones-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1]]])
  =/  canon-gth-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint tail=~] baum=~[~[0 0 0]]])
  =/  assay-gth-1x3-6u  (gth:la input-ones-1x3-6u jnput-ones-1x3-6u)
  %+  is-equal
    canon-gth-1x3-6u
  assay-gth-1x3-6u
::
::
++  test-gth-2x1-6u  ^-  tang
  =/  input-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  jnput-ones-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1]]])
  =/  canon-gth-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint tail=~] baum=~[~[0] ~[0]]])
  =/  assay-gth-2x1-6u  (gth:la input-ones-2x1-6u jnput-ones-2x1-6u)
  %+  is-equal
    canon-gth-2x1-6u
  assay-gth-2x1-6u
::
::
++  test-gth-2x2-6u  ^-  tang
  =/  input-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  jnput-ones-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1]]])
  =/  canon-gth-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint tail=~] baum=~[~[0 0] ~[0 0]]])
  =/  assay-gth-2x2-6u  (gth:la input-ones-2x2-6u jnput-ones-2x2-6u)
  %+  is-equal
    canon-gth-2x2-6u
  assay-gth-2x2-6u
::
::
++  test-gth-2x3-6u  ^-  tang
  =/  input-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1]]])
  =/  canon-gth-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0]]])
  =/  assay-gth-2x3-6u  (gth:la input-ones-2x3-6u jnput-ones-2x3-6u)
  %+  is-equal
    canon-gth-2x3-6u
  assay-gth-2x3-6u
::
::
++  test-gth-3x1-6u  ^-  tang
  =/  input-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  jnput-ones-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[1] ~[1] ~[1]]])
  =/  canon-gth-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint tail=~] baum=~[~[0] ~[0] ~[0]]])
  =/  assay-gth-3x1-6u  (gth:la input-ones-3x1-6u jnput-ones-3x1-6u)
  %+  is-equal
    canon-gth-3x1-6u
  assay-gth-3x1-6u
::
::
++  test-gth-3x2-6u  ^-  tang
  =/  input-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  jnput-ones-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[1 1] ~[1 1] ~[1 1]]])
  =/  canon-gth-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint tail=~] baum=~[~[0 0] ~[0 0] ~[0 0]]])
  =/  assay-gth-3x2-6u  (gth:la input-ones-3x2-6u jnput-ones-3x2-6u)
  %+  is-equal
    canon-gth-3x2-6u
  assay-gth-3x2-6u
::
::
++  test-gth-3x3-6u  ^-  tang
  =/  input-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  jnput-ones-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[1 1 1] ~[1 1 1] ~[1 1 1]]])
  =/  canon-gth-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint tail=~] baum=~[~[0 0 0] ~[0 0 0] ~[0 0 0]]])
  =/  assay-gth-3x3-6u  (gth:la input-ones-3x3-6u jnput-ones-3x3-6u)
  %+  is-equal
    canon-gth-3x3-6u
  assay-gth-3x3-6u
::

++  test-max-2-4r  ^-  tang
  =/  input-max-2-4r  (en-ray:la [meta=[shape=~[2] bloq=4 kind=%i754 prec=~] baum=(reap 2 .~~0.0)])
  =/  canon-max-2-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0]]])
  =/  assay-max-2-4r  (max:la (reshape:la (linspace:la meta.input-max-2-4r [.~~0.0 .~~1.0] 2) ~[1 2]))
  %+  is-equal
    canon-max-2-4r
  assay-max-2-4r
::

++  test-max-9-4r  ^-  tang
  =/  input-max-9-4r  (en-ray:la [meta=[shape=~[9] bloq=4 kind=%i754 prec=~] baum=(reap 9 .~~0.0)])
  =/  canon-max-9-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~1.0]]])
  =/  assay-max-9-4r  (max:la (reshape:la (linspace:la meta.input-max-9-4r [.~~0.0 .~~1.0] 9) ~[1 9]))
  %+  is-equal
    canon-max-9-4r
  assay-max-9-4r
::

++  test-max-2-5r  ^-  tang
  =/  input-max-2-5r  (en-ray:la [meta=[shape=~[2] bloq=5 kind=%i754 prec=~] baum=(reap 2 .0.0)])
  =/  canon-max-2-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0]]])
  =/  assay-max-2-5r  (max:la (reshape:la (linspace:la meta.input-max-2-5r [.0.0 .1.0] 2) ~[1 2]))
  %+  is-equal
    canon-max-2-5r
  assay-max-2-5r
::

++  test-max-9-5r  ^-  tang
  =/  input-max-9-5r  (en-ray:la [meta=[shape=~[9] bloq=5 kind=%i754 prec=~] baum=(reap 9 .0.0)])
  =/  canon-max-9-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 prec=~] baum=~[~[.1.0]]])
  =/  assay-max-9-5r  (max:la (reshape:la (linspace:la meta.input-max-9-5r [.0.0 .1.0] 9) ~[1 9]))
  %+  is-equal
    canon-max-9-5r
  assay-max-9-5r
::

++  test-max-2-6r  ^-  tang
  =/  input-max-2-6r  (en-ray:la [meta=[shape=~[2] bloq=6 kind=%i754 prec=~] baum=(reap 2 .~0.0)])
  =/  canon-max-2-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0]]])
  =/  assay-max-2-6r  (max:la (reshape:la (linspace:la meta.input-max-2-6r [.~0.0 .~1.0] 2) ~[1 2]))
  %+  is-equal
    canon-max-2-6r
  assay-max-2-6r
::

++  test-max-9-6r  ^-  tang
  =/  input-max-9-6r  (en-ray:la [meta=[shape=~[9] bloq=6 kind=%i754 prec=~] baum=(reap 9 .~0.0)])
  =/  canon-max-9-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~1.0]]])
  =/  assay-max-9-6r  (max:la (reshape:la (linspace:la meta.input-max-9-6r [.~0.0 .~1.0] 9) ~[1 9]))
  %+  is-equal
    canon-max-9-6r
  assay-max-9-6r
::

++  test-max-2-7r  ^-  tang
  =/  input-max-2-7r  (en-ray:la [meta=[shape=~[2] bloq=7 kind=%i754 prec=~] baum=(reap 2 .~~~0.0)])
  =/  canon-max-2-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0]]])
  =/  assay-max-2-7r  (max:la (reshape:la (linspace:la meta.input-max-2-7r [.~~~0.0 .~~~1.0] 2) ~[1 2]))
  %+  is-equal
    canon-max-2-7r
  assay-max-2-7r
::

++  test-max-9-7r  ^-  tang
  =/  input-max-9-7r  (en-ray:la [meta=[shape=~[9] bloq=7 kind=%i754 prec=~] baum=(reap 9 .~~~0.0)])
  =/  canon-max-9-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~1.0]]])
  =/  assay-max-9-7r  (max:la (reshape:la (linspace:la meta.input-max-9-7r [.~~~0.0 .~~~1.0] 9) ~[1 9]))
  %+  is-equal
    canon-max-9-7r
  assay-max-9-7r
::

++  test-max-1x1-3u  ^-  tang
  =/  input-max-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  canon-max-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-max-1x1-3u  (max:la (magic:la meta.input-max-1x1-3u))
  %+  is-equal
    canon-max-1x1-3u
  assay-max-1x1-3u
::

++  test-max-1x2-3u  ^-  tang
  =/  input-max-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-max-1x2-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[1]])
  =/  assay-max-1x2-3u  (max:la (magic:la meta.input-max-1x2-3u))
  %+  is-equal
    canon-max-1x2-3u
  assay-max-1x2-3u
::

++  test-max-1x3-3u  ^-  tang
  =/  input-max-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint prec=~] baum=~[~[0 1 2]]])
  =/  canon-max-1x3-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[2]])
  =/  assay-max-1x3-3u  (max:la (magic:la meta.input-max-1x3-3u))
  %+  is-equal
    canon-max-1x3-3u
  assay-max-1x3-3u
::

++  test-max-2x1-3u  ^-  tang
  =/  input-max-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  canon-max-2x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[1]])
  =/  assay-max-2x1-3u  (max:la (magic:la meta.input-max-2x1-3u))
  %+  is-equal
    canon-max-2x1-3u
  assay-max-2x1-3u
::

++  test-max-2x2-3u  ^-  tang
  =/  input-max-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint prec=~] baum=~[~[0 1] ~[2 3]]])
  =/  canon-max-2x2-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[3]])
  =/  assay-max-2x2-3u  (max:la (magic:la meta.input-max-2x2-3u))
  %+  is-equal
    canon-max-2x2-3u
  assay-max-2x2-3u
::

++  test-max-2x3-3u  ^-  tang
  =/  input-max-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5]]])
  =/  canon-max-2x3-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[5]])
  =/  assay-max-2x3-3u  (max:la (magic:la meta.input-max-2x3-3u))
  %+  is-equal
    canon-max-2x3-3u
  assay-max-2x3-3u
::

++  test-max-3x1-3u  ^-  tang
  =/  input-max-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint prec=~] baum=~[~[0] ~[0] ~[2]]])
  =/  canon-max-3x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[2]])
  =/  assay-max-3x1-3u  (max:la (magic:la meta.input-max-3x1-3u))
  %+  is-equal
    canon-max-3x1-3u
  assay-max-3x1-3u
::

++  test-max-3x2-3u  ^-  tang
  =/  input-max-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint prec=~] baum=~[~[0 1] ~[2 3] ~[4 5]]])
  =/  canon-max-3x2-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[5]])
  =/  assay-max-3x2-3u  (max:la (magic:la meta.input-max-3x2-3u))
  %+  is-equal
    canon-max-3x2-3u
  assay-max-3x2-3u
::

++  test-max-3x3-3u  ^-  tang
  =/  input-max-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5] ~[6 7 8]]])
  =/  canon-max-3x3-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[8]])
  =/  assay-max-3x3-3u  (max:la (magic:la meta.input-max-3x3-3u))
  %+  is-equal
    canon-max-3x3-3u
  assay-max-3x3-3u
::

++  test-max-1x1-4u  ^-  tang
  =/  input-max-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  canon-max-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-max-1x1-4u  (max:la (magic:la meta.input-max-1x1-4u))
  %+  is-equal
    canon-max-1x1-4u
  assay-max-1x1-4u
::

++  test-max-1x2-4u  ^-  tang
  =/  input-max-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-max-1x2-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[1]])
  =/  assay-max-1x2-4u  (max:la (magic:la meta.input-max-1x2-4u))
  %+  is-equal
    canon-max-1x2-4u
  assay-max-1x2-4u
::

++  test-max-1x3-4u  ^-  tang
  =/  input-max-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint prec=~] baum=~[~[0 1 2]]])
  =/  canon-max-1x3-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[2]])
  =/  assay-max-1x3-4u  (max:la (magic:la meta.input-max-1x3-4u))
  %+  is-equal
    canon-max-1x3-4u
  assay-max-1x3-4u
::

++  test-max-2x1-4u  ^-  tang
  =/  input-max-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  canon-max-2x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[1]])
  =/  assay-max-2x1-4u  (max:la (magic:la meta.input-max-2x1-4u))
  %+  is-equal
    canon-max-2x1-4u
  assay-max-2x1-4u
::

++  test-max-2x2-4u  ^-  tang
  =/  input-max-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint prec=~] baum=~[~[0 1] ~[2 3]]])
  =/  canon-max-2x2-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[3]])
  =/  assay-max-2x2-4u  (max:la (magic:la meta.input-max-2x2-4u))
  %+  is-equal
    canon-max-2x2-4u
  assay-max-2x2-4u
::

++  test-max-2x3-4u  ^-  tang
  =/  input-max-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5]]])
  =/  canon-max-2x3-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[5]])
  =/  assay-max-2x3-4u  (max:la (magic:la meta.input-max-2x3-4u))
  %+  is-equal
    canon-max-2x3-4u
  assay-max-2x3-4u
::

++  test-max-3x1-4u  ^-  tang
  =/  input-max-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint prec=~] baum=~[~[0] ~[0] ~[2]]])
  =/  canon-max-3x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[2]])
  =/  assay-max-3x1-4u  (max:la (magic:la meta.input-max-3x1-4u))
  %+  is-equal
    canon-max-3x1-4u
  assay-max-3x1-4u
::

++  test-max-3x2-4u  ^-  tang
  =/  input-max-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint prec=~] baum=~[~[0 1] ~[2 3] ~[4 5]]])
  =/  canon-max-3x2-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[5]])
  =/  assay-max-3x2-4u  (max:la (magic:la meta.input-max-3x2-4u))
  %+  is-equal
    canon-max-3x2-4u
  assay-max-3x2-4u
::

++  test-max-3x3-4u  ^-  tang
  =/  input-max-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5] ~[6 7 8]]])
  =/  canon-max-3x3-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[8]])
  =/  assay-max-3x3-4u  (max:la (magic:la meta.input-max-3x3-4u))
  %+  is-equal
    canon-max-3x3-4u
  assay-max-3x3-4u
::

++  test-max-1x1-5u  ^-  tang
  =/  input-max-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  canon-max-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-max-1x1-5u  (max:la (magic:la meta.input-max-1x1-5u))
  %+  is-equal
    canon-max-1x1-5u
  assay-max-1x1-5u
::

++  test-max-1x2-5u  ^-  tang
  =/  input-max-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-max-1x2-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[1]])
  =/  assay-max-1x2-5u  (max:la (magic:la meta.input-max-1x2-5u))
  %+  is-equal
    canon-max-1x2-5u
  assay-max-1x2-5u
::

++  test-max-1x3-5u  ^-  tang
  =/  input-max-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint prec=~] baum=~[~[0 1 2]]])
  =/  canon-max-1x3-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[2]])
  =/  assay-max-1x3-5u  (max:la (magic:la meta.input-max-1x3-5u))
  %+  is-equal
    canon-max-1x3-5u
  assay-max-1x3-5u
::

++  test-max-2x1-5u  ^-  tang
  =/  input-max-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  canon-max-2x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[1]])
  =/  assay-max-2x1-5u  (max:la (magic:la meta.input-max-2x1-5u))
  %+  is-equal
    canon-max-2x1-5u
  assay-max-2x1-5u
::

++  test-max-2x2-5u  ^-  tang
  =/  input-max-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint prec=~] baum=~[~[0 1] ~[2 3]]])
  =/  canon-max-2x2-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[3]])
  =/  assay-max-2x2-5u  (max:la (magic:la meta.input-max-2x2-5u))
  %+  is-equal
    canon-max-2x2-5u
  assay-max-2x2-5u
::

++  test-max-2x3-5u  ^-  tang
  =/  input-max-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5]]])
  =/  canon-max-2x3-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[5]])
  =/  assay-max-2x3-5u  (max:la (magic:la meta.input-max-2x3-5u))
  %+  is-equal
    canon-max-2x3-5u
  assay-max-2x3-5u
::

++  test-max-3x1-5u  ^-  tang
  =/  input-max-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint prec=~] baum=~[~[0] ~[0] ~[2]]])
  =/  canon-max-3x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[2]])
  =/  assay-max-3x1-5u  (max:la (magic:la meta.input-max-3x1-5u))
  %+  is-equal
    canon-max-3x1-5u
  assay-max-3x1-5u
::

++  test-max-3x2-5u  ^-  tang
  =/  input-max-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint prec=~] baum=~[~[0 1] ~[2 3] ~[4 5]]])
  =/  canon-max-3x2-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[5]])
  =/  assay-max-3x2-5u  (max:la (magic:la meta.input-max-3x2-5u))
  %+  is-equal
    canon-max-3x2-5u
  assay-max-3x2-5u
::

++  test-max-3x3-5u  ^-  tang
  =/  input-max-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5] ~[6 7 8]]])
  =/  canon-max-3x3-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[8]])
  =/  assay-max-3x3-5u  (max:la (magic:la meta.input-max-3x3-5u))
  %+  is-equal
    canon-max-3x3-5u
  assay-max-3x3-5u
::

++  test-max-1x1-6u  ^-  tang
  =/  input-max-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  canon-max-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-max-1x1-6u  (max:la (magic:la meta.input-max-1x1-6u))
  %+  is-equal
    canon-max-1x1-6u
  assay-max-1x1-6u
::

++  test-max-1x2-6u  ^-  tang
  =/  input-max-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-max-1x2-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[1]])
  =/  assay-max-1x2-6u  (max:la (magic:la meta.input-max-1x2-6u))
  %+  is-equal
    canon-max-1x2-6u
  assay-max-1x2-6u
::

++  test-max-1x3-6u  ^-  tang
  =/  input-max-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint prec=~] baum=~[~[0 1 2]]])
  =/  canon-max-1x3-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[2]])
  =/  assay-max-1x3-6u  (max:la (magic:la meta.input-max-1x3-6u))
  %+  is-equal
    canon-max-1x3-6u
  assay-max-1x3-6u
::

++  test-max-2x1-6u  ^-  tang
  =/  input-max-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  canon-max-2x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[1]])
  =/  assay-max-2x1-6u  (max:la (magic:la meta.input-max-2x1-6u))
  %+  is-equal
    canon-max-2x1-6u
  assay-max-2x1-6u
::

++  test-max-2x2-6u  ^-  tang
  =/  input-max-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint prec=~] baum=~[~[0 1] ~[2 3]]])
  =/  canon-max-2x2-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[3]])
  =/  assay-max-2x2-6u  (max:la (magic:la meta.input-max-2x2-6u))
  %+  is-equal
    canon-max-2x2-6u
  assay-max-2x2-6u
::

++  test-max-2x3-6u  ^-  tang
  =/  input-max-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5]]])
  =/  canon-max-2x3-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[5]])
  =/  assay-max-2x3-6u  (max:la (magic:la meta.input-max-2x3-6u))
  %+  is-equal
    canon-max-2x3-6u
  assay-max-2x3-6u
::

++  test-max-3x1-6u  ^-  tang
  =/  input-max-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint prec=~] baum=~[~[0] ~[0] ~[2]]])
  =/  canon-max-3x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[2]])
  =/  assay-max-3x1-6u  (max:la (magic:la meta.input-max-3x1-6u))
  %+  is-equal
    canon-max-3x1-6u
  assay-max-3x1-6u
::

++  test-max-3x2-6u  ^-  tang
  =/  input-max-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint prec=~] baum=~[~[0 1] ~[2 3] ~[4 5]]])
  =/  canon-max-3x2-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[5]])
  =/  assay-max-3x2-6u  (max:la (magic:la meta.input-max-3x2-6u))
  %+  is-equal
    canon-max-3x2-6u
  assay-max-3x2-6u
::

++  test-max-3x3-6u  ^-  tang
  =/  input-max-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5] ~[6 7 8]]])
  =/  canon-max-3x3-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[8]])
  =/  assay-max-3x3-6u  (max:la (magic:la meta.input-max-3x3-6u))
  %+  is-equal
    canon-max-3x3-6u
  assay-max-3x3-6u
::

++  test-min-2-4r  ^-  tang
  =/  input-min-2-4r  (en-ray:la [meta=[shape=~[2] bloq=4 kind=%i754 prec=~] baum=(reap 2 .~~0.0)])
  =/  canon-min-2-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0]]])
  =/  assay-min-2-4r  (min:la (reshape:la (linspace:la meta.input-min-2-4r [.~~0.0 .~~1.0] 2) ~[1 2]))
  %+  is-equal
    canon-min-2-4r
  assay-min-2-4r
::

++  test-min-9-4r  ^-  tang
  =/  input-min-9-4r  (en-ray:la [meta=[shape=~[9] bloq=4 kind=%i754 prec=~] baum=(reap 9 .~~0.0)])
  =/  canon-min-9-4r  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%i754 prec=~] baum=~[~[.~~0.0]]])
  =/  assay-min-9-4r  (min:la (reshape:la (linspace:la meta.input-min-9-4r [.~~0.0 .~~1.0] 9) ~[1 9]))
  %+  is-equal
    canon-min-9-4r
  assay-min-9-4r
::

++  test-min-2-5r  ^-  tang
  =/  input-min-2-5r  (en-ray:la [meta=[shape=~[2] bloq=5 kind=%i754 prec=~] baum=(reap 2 .0.0)])
  =/  canon-min-2-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0]]])
  =/  assay-min-2-5r  (min:la (reshape:la (linspace:la meta.input-min-2-5r [.0.0 .1.0] 2) ~[1 2]))
  %+  is-equal
    canon-min-2-5r
  assay-min-2-5r
::

++  test-min-9-5r  ^-  tang
  =/  input-min-9-5r  (en-ray:la [meta=[shape=~[9] bloq=5 kind=%i754 prec=~] baum=(reap 9 .0.0)])
  =/  canon-min-9-5r  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%i754 prec=~] baum=~[~[.0.0]]])
  =/  assay-min-9-5r  (min:la (reshape:la (linspace:la meta.input-min-9-5r [.0.0 .1.0] 9) ~[1 9]))
  %+  is-equal
    canon-min-9-5r
  assay-min-9-5r
::

++  test-min-2-6r  ^-  tang
  =/  input-min-2-6r  (en-ray:la [meta=[shape=~[2] bloq=6 kind=%i754 prec=~] baum=(reap 2 .~0.0)])
  =/  canon-min-2-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0]]])
  =/  assay-min-2-6r  (min:la (reshape:la (linspace:la meta.input-min-2-6r [.~0.0 .~1.0] 2) ~[1 2]))
  %+  is-equal
    canon-min-2-6r
  assay-min-2-6r
::

++  test-min-9-6r  ^-  tang
  =/  input-min-9-6r  (en-ray:la [meta=[shape=~[9] bloq=6 kind=%i754 prec=~] baum=(reap 9 .~0.0)])
  =/  canon-min-9-6r  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%i754 prec=~] baum=~[~[.~0.0]]])
  =/  assay-min-9-6r  (min:la (reshape:la (linspace:la meta.input-min-9-6r [.~0.0 .~1.0] 9) ~[1 9]))
  %+  is-equal
    canon-min-9-6r
  assay-min-9-6r
::

++  test-min-2-7r  ^-  tang
  =/  input-min-2-7r  (en-ray:la [meta=[shape=~[2] bloq=7 kind=%i754 prec=~] baum=(reap 2 .~~~0.0)])
  =/  canon-min-2-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0]]])
  =/  assay-min-2-7r  (min:la (reshape:la (linspace:la meta.input-min-2-7r [.~~~0.0 .~~~1.0] 2) ~[1 2]))
  %+  is-equal
    canon-min-2-7r
  assay-min-2-7r
::

++  test-min-9-7r  ^-  tang
  =/  input-min-9-7r  (en-ray:la [meta=[shape=~[9] bloq=7 kind=%i754 prec=~] baum=(reap 9 .~~~0.0)])
  =/  canon-min-9-7r  (en-ray:la [meta=[shape=~[1 1] bloq=7 kind=%i754 prec=~] baum=~[~[.~~~0.0]]])
  =/  assay-min-9-7r  (min:la (reshape:la (linspace:la meta.input-min-9-7r [.~~~0.0 .~~~1.0] 9) ~[1 9]))
  %+  is-equal
    canon-min-9-7r
  assay-min-9-7r
::

++  test-min-1x1-3u  ^-  tang
  =/  input-min-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  canon-min-1x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x1-3u  (min:la (magic:la meta.input-min-1x1-3u))
  %+  is-equal
    canon-min-1x1-3u
  assay-min-1x1-3u
::

++  test-min-1x2-3u  ^-  tang
  =/  input-min-1x2-3u  (en-ray:la [meta=[shape=~[1 2] bloq=3 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-min-1x2-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x2-3u  (min:la (magic:la meta.input-min-1x2-3u))
  %+  is-equal
    canon-min-1x2-3u
  assay-min-1x2-3u
::

++  test-min-1x3-3u  ^-  tang
  =/  input-min-1x3-3u  (en-ray:la [meta=[shape=~[1 3] bloq=3 kind=%uint prec=~] baum=~[~[0 1 2]]])
  =/  canon-min-1x3-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x3-3u  (min:la (magic:la meta.input-min-1x3-3u))
  %+  is-equal
    canon-min-1x3-3u
  assay-min-1x3-3u
::

++  test-min-2x1-3u  ^-  tang
  =/  input-min-2x1-3u  (en-ray:la [meta=[shape=~[2 1] bloq=3 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  canon-min-2x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x1-3u  (min:la (magic:la meta.input-min-2x1-3u))
  %+  is-equal
    canon-min-2x1-3u
  assay-min-2x1-3u
::

++  test-min-2x2-3u  ^-  tang
  =/  input-min-2x2-3u  (en-ray:la [meta=[shape=~[2 2] bloq=3 kind=%uint prec=~] baum=~[~[0 1] ~[2 3]]])
  =/  canon-min-2x2-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x2-3u  (min:la (magic:la meta.input-min-2x2-3u))
  %+  is-equal
    canon-min-2x2-3u
  assay-min-2x2-3u
::

++  test-min-2x3-3u  ^-  tang
  =/  input-min-2x3-3u  (en-ray:la [meta=[shape=~[2 3] bloq=3 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5]]])
  =/  canon-min-2x3-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x3-3u  (min:la (magic:la meta.input-min-2x3-3u))
  %+  is-equal
    canon-min-2x3-3u
  assay-min-2x3-3u
::

++  test-min-3x1-3u  ^-  tang
  =/  input-min-3x1-3u  (en-ray:la [meta=[shape=~[3 1] bloq=3 kind=%uint prec=~] baum=~[~[0] ~[0] ~[2]]])
  =/  canon-min-3x1-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x1-3u  (min:la (magic:la meta.input-min-3x1-3u))
  %+  is-equal
    canon-min-3x1-3u
  assay-min-3x1-3u
::

++  test-min-3x2-3u  ^-  tang
  =/  input-min-3x2-3u  (en-ray:la [meta=[shape=~[3 2] bloq=3 kind=%uint prec=~] baum=~[~[0 1] ~[2 3] ~[4 5]]])
  =/  canon-min-3x2-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x2-3u  (min:la (magic:la meta.input-min-3x2-3u))
  %+  is-equal
    canon-min-3x2-3u
  assay-min-3x2-3u
::

++  test-min-3x3-3u  ^-  tang
  =/  input-min-3x3-3u  (en-ray:la [meta=[shape=~[3 3] bloq=3 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5] ~[6 7 8]]])
  =/  canon-min-3x3-3u  (en-ray:la [meta=[shape=~[1 1] bloq=3 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x3-3u  (min:la (magic:la meta.input-min-3x3-3u))
  %+  is-equal
    canon-min-3x3-3u
  assay-min-3x3-3u
::

++  test-min-1x1-4u  ^-  tang
  =/  input-min-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  canon-min-1x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x1-4u  (min:la (magic:la meta.input-min-1x1-4u))
  %+  is-equal
    canon-min-1x1-4u
  assay-min-1x1-4u
::

++  test-min-1x2-4u  ^-  tang
  =/  input-min-1x2-4u  (en-ray:la [meta=[shape=~[1 2] bloq=4 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-min-1x2-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x2-4u  (min:la (magic:la meta.input-min-1x2-4u))
  %+  is-equal
    canon-min-1x2-4u
  assay-min-1x2-4u
::

++  test-min-1x3-4u  ^-  tang
  =/  input-min-1x3-4u  (en-ray:la [meta=[shape=~[1 3] bloq=4 kind=%uint prec=~] baum=~[~[0 1 2]]])
  =/  canon-min-1x3-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x3-4u  (min:la (magic:la meta.input-min-1x3-4u))
  %+  is-equal
    canon-min-1x3-4u
  assay-min-1x3-4u
::

++  test-min-2x1-4u  ^-  tang
  =/  input-min-2x1-4u  (en-ray:la [meta=[shape=~[2 1] bloq=4 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  canon-min-2x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x1-4u  (min:la (magic:la meta.input-min-2x1-4u))
  %+  is-equal
    canon-min-2x1-4u
  assay-min-2x1-4u
::

++  test-min-2x2-4u  ^-  tang
  =/  input-min-2x2-4u  (en-ray:la [meta=[shape=~[2 2] bloq=4 kind=%uint prec=~] baum=~[~[0 1] ~[2 3]]])
  =/  canon-min-2x2-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x2-4u  (min:la (magic:la meta.input-min-2x2-4u))
  %+  is-equal
    canon-min-2x2-4u
  assay-min-2x2-4u
::

++  test-min-2x3-4u  ^-  tang
  =/  input-min-2x3-4u  (en-ray:la [meta=[shape=~[2 3] bloq=4 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5]]])
  =/  canon-min-2x3-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x3-4u  (min:la (magic:la meta.input-min-2x3-4u))
  %+  is-equal
    canon-min-2x3-4u
  assay-min-2x3-4u
::

++  test-min-3x1-4u  ^-  tang
  =/  input-min-3x1-4u  (en-ray:la [meta=[shape=~[3 1] bloq=4 kind=%uint prec=~] baum=~[~[0] ~[0] ~[2]]])
  =/  canon-min-3x1-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x1-4u  (min:la (magic:la meta.input-min-3x1-4u))
  %+  is-equal
    canon-min-3x1-4u
  assay-min-3x1-4u
::

++  test-min-3x2-4u  ^-  tang
  =/  input-min-3x2-4u  (en-ray:la [meta=[shape=~[3 2] bloq=4 kind=%uint prec=~] baum=~[~[0 1] ~[2 3] ~[4 5]]])
  =/  canon-min-3x2-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x2-4u  (min:la (magic:la meta.input-min-3x2-4u))
  %+  is-equal
    canon-min-3x2-4u
  assay-min-3x2-4u
::

++  test-min-3x3-4u  ^-  tang
  =/  input-min-3x3-4u  (en-ray:la [meta=[shape=~[3 3] bloq=4 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5] ~[6 7 8]]])
  =/  canon-min-3x3-4u  (en-ray:la [meta=[shape=~[1 1] bloq=4 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x3-4u  (min:la (magic:la meta.input-min-3x3-4u))
  %+  is-equal
    canon-min-3x3-4u
  assay-min-3x3-4u
::

++  test-min-1x1-5u  ^-  tang
  =/  input-min-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  canon-min-1x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x1-5u  (min:la (magic:la meta.input-min-1x1-5u))
  %+  is-equal
    canon-min-1x1-5u
  assay-min-1x1-5u
::

++  test-min-1x2-5u  ^-  tang
  =/  input-min-1x2-5u  (en-ray:la [meta=[shape=~[1 2] bloq=5 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-min-1x2-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x2-5u  (min:la (magic:la meta.input-min-1x2-5u))
  %+  is-equal
    canon-min-1x2-5u
  assay-min-1x2-5u
::

++  test-min-1x3-5u  ^-  tang
  =/  input-min-1x3-5u  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%uint prec=~] baum=~[~[0 1 2]]])
  =/  canon-min-1x3-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x3-5u  (min:la (magic:la meta.input-min-1x3-5u))
  %+  is-equal
    canon-min-1x3-5u
  assay-min-1x3-5u
::

++  test-min-2x1-5u  ^-  tang
  =/  input-min-2x1-5u  (en-ray:la [meta=[shape=~[2 1] bloq=5 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  canon-min-2x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x1-5u  (min:la (magic:la meta.input-min-2x1-5u))
  %+  is-equal
    canon-min-2x1-5u
  assay-min-2x1-5u
::

++  test-min-2x2-5u  ^-  tang
  =/  input-min-2x2-5u  (en-ray:la [meta=[shape=~[2 2] bloq=5 kind=%uint prec=~] baum=~[~[0 1] ~[2 3]]])
  =/  canon-min-2x2-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x2-5u  (min:la (magic:la meta.input-min-2x2-5u))
  %+  is-equal
    canon-min-2x2-5u
  assay-min-2x2-5u
::

++  test-min-2x3-5u  ^-  tang
  =/  input-min-2x3-5u  (en-ray:la [meta=[shape=~[2 3] bloq=5 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5]]])
  =/  canon-min-2x3-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x3-5u  (min:la (magic:la meta.input-min-2x3-5u))
  %+  is-equal
    canon-min-2x3-5u
  assay-min-2x3-5u
::

++  test-min-3x1-5u  ^-  tang
  =/  input-min-3x1-5u  (en-ray:la [meta=[shape=~[3 1] bloq=5 kind=%uint prec=~] baum=~[~[0] ~[0] ~[2]]])
  =/  canon-min-3x1-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x1-5u  (min:la (magic:la meta.input-min-3x1-5u))
  %+  is-equal
    canon-min-3x1-5u
  assay-min-3x1-5u
::

++  test-min-3x2-5u  ^-  tang
  =/  input-min-3x2-5u  (en-ray:la [meta=[shape=~[3 2] bloq=5 kind=%uint prec=~] baum=~[~[0 1] ~[2 3] ~[4 5]]])
  =/  canon-min-3x2-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x2-5u  (min:la (magic:la meta.input-min-3x2-5u))
  %+  is-equal
    canon-min-3x2-5u
  assay-min-3x2-5u
::

++  test-min-3x3-5u  ^-  tang
  =/  input-min-3x3-5u  (en-ray:la [meta=[shape=~[3 3] bloq=5 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5] ~[6 7 8]]])
  =/  canon-min-3x3-5u  (en-ray:la [meta=[shape=~[1 1] bloq=5 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x3-5u  (min:la (magic:la meta.input-min-3x3-5u))
  %+  is-equal
    canon-min-3x3-5u
  assay-min-3x3-5u
::

++  test-min-1x1-6u  ^-  tang
  =/  input-min-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  canon-min-1x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x1-6u  (min:la (magic:la meta.input-min-1x1-6u))
  %+  is-equal
    canon-min-1x1-6u
  assay-min-1x1-6u
::

++  test-min-1x2-6u  ^-  tang
  =/  input-min-1x2-6u  (en-ray:la [meta=[shape=~[1 2] bloq=6 kind=%uint prec=~] baum=~[~[1 1]]])
  =/  canon-min-1x2-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x2-6u  (min:la (magic:la meta.input-min-1x2-6u))
  %+  is-equal
    canon-min-1x2-6u
  assay-min-1x2-6u
::

++  test-min-1x3-6u  ^-  tang
  =/  input-min-1x3-6u  (en-ray:la [meta=[shape=~[1 3] bloq=6 kind=%uint prec=~] baum=~[~[0 1 2]]])
  =/  canon-min-1x3-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-1x3-6u  (min:la (magic:la meta.input-min-1x3-6u))
  %+  is-equal
    canon-min-1x3-6u
  assay-min-1x3-6u
::

++  test-min-2x1-6u  ^-  tang
  =/  input-min-2x1-6u  (en-ray:la [meta=[shape=~[2 1] bloq=6 kind=%uint prec=~] baum=~[~[0] ~[0]]])
  =/  canon-min-2x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x1-6u  (min:la (magic:la meta.input-min-2x1-6u))
  %+  is-equal
    canon-min-2x1-6u
  assay-min-2x1-6u
::

++  test-min-2x2-6u  ^-  tang
  =/  input-min-2x2-6u  (en-ray:la [meta=[shape=~[2 2] bloq=6 kind=%uint prec=~] baum=~[~[0 1] ~[2 3]]])
  =/  canon-min-2x2-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x2-6u  (min:la (magic:la meta.input-min-2x2-6u))
  %+  is-equal
    canon-min-2x2-6u
  assay-min-2x2-6u
::

++  test-min-2x3-6u  ^-  tang
  =/  input-min-2x3-6u  (en-ray:la [meta=[shape=~[2 3] bloq=6 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5]]])
  =/  canon-min-2x3-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-2x3-6u  (min:la (magic:la meta.input-min-2x3-6u))
  %+  is-equal
    canon-min-2x3-6u
  assay-min-2x3-6u
::

++  test-min-3x1-6u  ^-  tang
  =/  input-min-3x1-6u  (en-ray:la [meta=[shape=~[3 1] bloq=6 kind=%uint prec=~] baum=~[~[0] ~[0] ~[2]]])
  =/  canon-min-3x1-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x1-6u  (min:la (magic:la meta.input-min-3x1-6u))
  %+  is-equal
    canon-min-3x1-6u
  assay-min-3x1-6u
::

++  test-min-3x2-6u  ^-  tang
  =/  input-min-3x2-6u  (en-ray:la [meta=[shape=~[3 2] bloq=6 kind=%uint prec=~] baum=~[~[0 1] ~[2 3] ~[4 5]]])
  =/  canon-min-3x2-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x2-6u  (min:la (magic:la meta.input-min-3x2-6u))
  %+  is-equal
    canon-min-3x2-6u
  assay-min-3x2-6u
::

++  test-min-3x3-6u  ^-  tang
  =/  input-min-3x3-6u  (en-ray:la [meta=[shape=~[3 3] bloq=6 kind=%uint prec=~] baum=~[~[0 1 2] ~[3 4 5] ~[6 7 8]]])
  =/  canon-min-3x3-6u  (en-ray:la [meta=[shape=~[1 1] bloq=6 kind=%uint prec=~] baum=~[~[0]]])
  =/  assay-min-3x3-6u  (min:la (magic:la meta.input-min-3x3-6u))
  %+  is-equal
    canon-min-3x3-6u
  assay-min-3x3-6u
::
::
::  argmin/argmax return the ravel index of the min/max element (matching
::  the Hoon `find` reference), so the element fetched at that index must
::  equal the min/max value.  A reversed index (len - i - 1) would fetch
::  the wrong element for a non-symmetric array.
++  test-argmin-argmax  ^-  tang
  =/  a  (en-ray:la [meta=[shape=~[1 4] bloq=5 kind=%i754 tail=~] baum=~[~[.3.0 .1.0 .4.0 .2.0]]])
  ;:  weld
    %+  expect-eq
      !>((get-item:la (min:la a) ~[0 0]))
      !>((snag (argmin:la a) (ravel:la a)))
    %+  expect-eq
      !>((get-item:la (max:la a) ~[0 0]))
      !>((snag (argmax:la a) (ravel:la a)))
  ==
::
::
::  any/all over boolean rays (numeric convention: true=1.0, false=0.0).
::  any = some element truthy; all = every element truthy.
++  test-any-all  ^-  tang
  =/  all-true   (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .1.0 .1.0]]])
  =/  has-false  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.1.0 .0.0 .1.0]]])
  =/  all-false  (en-ray:la [meta=[shape=~[1 3] bloq=5 kind=%i754 tail=~] baum=~[~[.0.0 .0.0 .0.0]]])
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((all:la all-true))
    %+  expect-eq  !>(%.y)  !>((any:la all-true))
    %+  expect-eq  !>(%.n)  !>((all:la has-false))
    %+  expect-eq  !>(%.y)  !>((any:la has-false))
    %+  expect-eq  !>(%.n)  !>((all:la all-false))
    %+  expect-eq  !>(%.n)  !>((any:la all-false))
  ==
::
::
::  1-D coverage: min/max reduce regardless of rank, so any/all must work
::  on a plain vector too (the 2-D-only test above hid a rank assumption).
++  test-any-all-1d  ^-  tang
  =/  all-true   (en-ray:la [meta=[shape=~[3] bloq=5 kind=%i754 tail=~] baum=~[.1.0 .1.0 .1.0]])
  =/  has-false  (en-ray:la [meta=[shape=~[3] bloq=5 kind=%i754 tail=~] baum=~[.1.0 .0.0 .1.0]])
  =/  all-false  (en-ray:la [meta=[shape=~[3] bloq=5 kind=%i754 tail=~] baum=~[.0.0 .0.0 .0.0]])
  ;:  weld
    %+  expect-eq  !>(%.y)  !>((all:la all-true))
    %+  expect-eq  !>(%.y)  !>((any:la all-true))
    %+  expect-eq  !>(%.n)  !>((all:la has-false))
    %+  expect-eq  !>(%.y)  !>((any:la has-false))
    %+  expect-eq  !>(%.n)  !>((all:la all-false))
    %+  expect-eq  !>(%.n)  !>((any:la all-false))
  ==
--
