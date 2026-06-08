::  Hermitian eig at half (@ch) and quad (@cq) width: [[2,i],[-i,2]] -> {1,3}.
/-  ls=lagoon
/+  *test, *saloon, *lagoon, complex
|%
++  sch  (sake %n .~~1e-2)
++  ach
  =/  pk  ~(pak ch:complex %n)
  =/  cj  ~(conj ch:complex %n)
  ^-  ray:ls
  (en-ray:(lake %n) [[~[2 2] 5 %cplx ~] ~[~[(pk .~~2 .~~0) (pk .~~0 .~~1)] ~[(cj (pk .~~0 .~~1)) (pk .~~2 .~~0)]]])
++  test-ch
  %-  expect
  !>  (all:(lake %n) (is-close:(lake %n) (eigvals:sch ach) (en-ray:(lake %n) [[~[2] 4 %i754 ~] ~[.~~1 .~~3]]) [.~~1e-1 .~~1e-1]))
++  scq  (sake %n .~~~1e-25)
++  acq
  =/  pk  ~(pak cq:complex %n)
  =/  cj  ~(conj cq:complex %n)
  ^-  ray:ls
  (en-ray:(lake %n) [[~[2 2] 8 %cplx ~] ~[~[(pk .~~~2 .~~~0) (pk .~~~0 .~~~1)] ~[(cj (pk .~~~0 .~~~1)) (pk .~~~2 .~~~0)]]])
++  test-cq
  %-  expect
  !>  (all:(lake %n) (is-close:(lake %n) (eigvals:scq acq) (en-ray:(lake %n) [[~[2] 7 %i754 ~] ~[.~~~1 .~~~3]]) [.~~~1e-20 .~~~1e-20]))
--
