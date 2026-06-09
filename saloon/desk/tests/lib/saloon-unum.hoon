::  Saloon Tier-B transcendentals over %unum (posit) rays.  Verifies the
::  trans-scalar/fun-scalar %unum dispatch reaches /lib/unum per width.  Expected
::  scalar values are the same series outputs checked in /tests/lib/unum-fns
::  (tools/posit_check.py replica; correctly rounded vs mpmath at posit8).
::  Inputs at rps (posit32, bloq 5): 0x3800.0000 = .5.
::
/-  ls=lagoon
/+  *test, *saloon, *lagoon
|%
::  a 1-element %unum ray at bloq .b holding posit pattern .p
++  ur  |=([b=@ p=@] ^-(ray:ls (en-ray:(lake %n) [[~[1] b %unum ~] ~[p]])))
::  head scalar of a ray
++  hd  |=(a=ray:ls ^-(@ -:(ravel:(lake %n) a)))
++  test-unum-exp-rps   (expect-eq !>(`@`0x4530.94c8) !>((hd (exp:sa (ur 5 0x3800.0000)))))
++  test-unum-sin-rps   (expect-eq !>(`@`0x3757.743a) !>((hd (sin:sa (ur 5 0x3800.0000)))))
++  test-unum-cos-rps   (expect-eq !>(`@`0x3e0a.9404) !>((hd (cos:sa (ur 5 0x3800.0000)))))
++  test-unum-tan-rps   (expect-eq !>(`@`0x38bd.a7ad) !>((hd (tan:sa (ur 5 0x3800.0000)))))
++  test-unum-log-rps   (expect-eq !>(`@`0x3b17.2180) !>((hd (log:sa (ur 5 (sun:rps:unum 2))))))
++  test-unum-cbrt-rps  (expect-eq !>(`@`0x4000.0000) !>((hd (cbrt:sa (ur 5 (sun:rps:unum 1))))))
++  test-unum-fact-rps  (expect-eq !>(`@`0x5400.0000) !>((hd (factorial:sa (ur 5 (sun:rps:unum 3))))))
::  +pow-n is element-wise a^b over two rays; the exponent ray holds the raw
::  integer (read as @u by +fun-scalar), so its slot is a bare 3, not posit 3.
++  test-unum-pown-rps  (expect-eq !>(`@`0x5800.0000) !>((hd (pow-n:sa (ur 5 (sun:rps:unum 2)) (ur 5 3)))))
--
