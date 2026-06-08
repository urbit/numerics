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
--
