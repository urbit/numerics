::  eig accepts matrices symmetric/Hermitian WITHIN tolerance, rejects genuinely
::  asymmetric ones.  Guards the +near/+cnear relaxation (Gram-from-mmul / +-0).
/-  ls=lagoon
/+  *test, *saloon, *lagoon
|%
++  sad  (sake %n .~1e-12)
::  [[2, 1], [1+1ULP, 2]] -- symmetric up to one ULP; accepted, eigenvalues {1,3}.
++  a-near
  ^-  ray:ls
  (en-ray:(lake %n) [[~[2 2] 6 %i754 ~] ~[~[.~2 0x3ff0.0000.0000.0000] ~[0x3ff0.0000.0000.0001 .~2]]])
++  test-near-symmetric-ok
  %-  expect
  !>  (all:(lake %n) (is-close:(lake %n) (eigvals:sad a-near) (en-ray:(lake %n) [[~[2] 6 %i754 ~] ~[.~1 .~3]]) [.~1e-6 .~1e-6]))
::  [[2, 1], [2, 2]] -- off-diagonals differ by O(1); rejected (crashes).
++  a-asym
  ^-  ray:ls
  (en-ray:(lake %n) [[~[2 2] 6 %i754 ~] ~[~[.~2 .~1] ~[.~2 .~2]]])
++  test-asymmetric-crashes  (expect-fail |.((eigvals:sad a-asym)))
::  rtol width must match the component.  .~1e-12 is @rd (8 bytes) on an @rs
::  (bloq 5, 4-byte) array: wider than the component, so its high bytes would
::  silently mis-scale the threshold -- rejected.  A matching @rs rtol (.1e-6)
::  is accepted (and yields the two eigenvalues).
++  a-rs
  ^-  ray:ls
  (en-ray:(lake %n) [[~[2 2] 5 %i754 ~] ~[~[.2 .1] ~[.1 .2]]])
++  test-rtol-width-mismatch-rejected  (expect-fail |.((eigvals:(sake %n .~1e-12) a-rs)))
++  test-rtol-width-match-ok
  %-  expect
  !>  =(2 (lent (ravel:(lake %n) (eigvals:(sake %n .1e-6) a-rs))))
::  bare `sa` (no +sake): the unusable 0x1 default rtol is replaced by a width-
::  appropriate +feps, so eig still converges -- here to the same {1,3} as a
::  +sake-set run, rather than spinning out the 60-sweep cap.
++  test-bare-eig-converges
  =/  a  (en-ray:(lake %n) [[~[2 2] 6 %i754 ~] ~[~[.~2 .~1] ~[.~1 .~2]]])
  %-  expect
  !>  (all:(lake %n) (is-close:(lake %n) (eigvals:sa a) (eigvals:sad a) [.~1e-6 .~1e-6]))
--
