=/  rhs  `(list @rh)`~[`@rh`0x0 `@rh`0x8000 `@rh`0x3c00 `@rh`0xbc00 `@rh`0x4200 `@rh`0x2e66 `@rh`0x7c00 `@rh`0xfc00 `@rh`0x7e00 `@rh`0x3800 `@rh`0x4000 `@rh`0x3c01 `@rh`0xc200 `@rh`0x1 `@rh`0x7bff `@rh`0x400 `@rh`0x7bfe `@rh`0xfbff `@rh`0x8001 `@rh`0x3ff]
=/  rss  `(list @rs)`~[`@rs`0x0 `@rs`0x8000.0000 `@rs`0x3f80.0000 `@rs`0xbf80.0000 `@rs`0x4040.0000 `@rs`0x3dcc.cccd `@rs`0x7f80.0000 `@rs`0xff80.0000 `@rs`0x7fc0.0000 `@rs`0x3f00.0000 `@rs`0x4000.0000 `@rs`0x3f80.0001 `@rs`0xc040.0000 `@rs`0x1 `@rs`0x7f7f.ffff `@rs`0x80.0000 `@rs`0x7f7f.fffe `@rs`0xff7f.ffff `@rs`0x8000.0001 `@rs`0x7f.ffff]
=/  rds  `(list @rd)`~[`@rd`0x0 `@rd`0x8000.0000.0000.0000 `@rd`0x3ff0.0000.0000.0000 `@rd`0xbff0.0000.0000.0000 `@rd`0x4008.0000.0000.0000 `@rd`0x3fb9.9999.9999.999a `@rd`0x7ff0.0000.0000.0000 `@rd`0xfff0.0000.0000.0000 `@rd`0x7ff8.0000.0000.0000 `@rd`0x3fe0.0000.0000.0000 `@rd`0x4000.0000.0000.0000 `@rd`0x3ff0.0000.0000.0001 `@rd`0xc008.0000.0000.0000 `@rd`0x1 `@rd`0x7fef.ffff.ffff.ffff `@rd`0x10.0000.0000.0000 `@rd`0x7fef.ffff.ffff.fffe `@rd`0xffef.ffff.ffff.ffff `@rd`0x8000.0000.0000.0001 `@rd`0xf.ffff.ffff.ffff]
=/  rqs  `(list @rq)`~[`@rq`0x0 `@rq`0x8000.0000.0000.0000.0000.0000.0000.0000 `@rq`0x3fff.0000.0000.0000.0000.0000.0000.0000 `@rq`0xbfff.0000.0000.0000.0000.0000.0000.0000 `@rq`0x4000.8000.0000.0000.0000.0000.0000.0000 `@rq`0x3ffb.9999.9999.9999.9999.9999.9999.999a `@rq`0x7fff.0000.0000.0000.0000.0000.0000.0000 `@rq`0xffff.0000.0000.0000.0000.0000.0000.0000 `@rq`0x7fff.8000.0000.0000.0000.0000.0000.0000 `@rq`0x3ffe.0000.0000.0000.0000.0000.0000.0000 `@rq`0x4000.0000.0000.0000.0000.0000.0000.0000 `@rq`0x3fff.0000.0000.0000.0000.0000.0000.0001 `@rq`0xc000.8000.0000.0000.0000.0000.0000.0000 `@rq`0x1 `@rq`0x7ffe.ffff.ffff.ffff.ffff.ffff.ffff.ffff `@rq`0x1.0000.0000.0000.0000.0000.0000.0000 `@rq`0x7ffe.ffff.ffff.ffff.ffff.ffff.ffff.fffe `@rq`0xfffe.ffff.ffff.ffff.ffff.ffff.ffff.ffff `@rq`0x8000.0000.0000.0000.0000.0000.0000.0001 `@rq`0xffff.ffff.ffff.ffff.ffff.ffff.ffff]
;:  weld
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rhs
  |=  a=@rh
  %+  turn  rhs
  |=  b=@rh
  [%rh %n a b (~(add rh %n) a b) (~(sub rh %n) a b) (~(mul rh %n) a b) (~(div rh %n) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rhs
  |=  a=@rh
  %+  turn  rhs
  |=  b=@rh
  [%rh %u a b (~(add rh %u) a b) (~(sub rh %u) a b) (~(mul rh %u) a b) (~(div rh %u) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rhs
  |=  a=@rh
  %+  turn  rhs
  |=  b=@rh
  [%rh %d a b (~(add rh %d) a b) (~(sub rh %d) a b) (~(mul rh %d) a b) (~(div rh %d) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rhs
  |=  a=@rh
  %+  turn  rhs
  |=  b=@rh
  [%rh %z a b (~(add rh %z) a b) (~(sub rh %z) a b) (~(mul rh %z) a b) (~(div rh %z) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rss
  |=  a=@rs
  %+  turn  rss
  |=  b=@rs
  [%rs %n a b (~(add rs %n) a b) (~(sub rs %n) a b) (~(mul rs %n) a b) (~(div rs %n) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rss
  |=  a=@rs
  %+  turn  rss
  |=  b=@rs
  [%rs %u a b (~(add rs %u) a b) (~(sub rs %u) a b) (~(mul rs %u) a b) (~(div rs %u) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rss
  |=  a=@rs
  %+  turn  rss
  |=  b=@rs
  [%rs %d a b (~(add rs %d) a b) (~(sub rs %d) a b) (~(mul rs %d) a b) (~(div rs %d) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rss
  |=  a=@rs
  %+  turn  rss
  |=  b=@rs
  [%rs %z a b (~(add rs %z) a b) (~(sub rs %z) a b) (~(mul rs %z) a b) (~(div rs %z) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rds
  |=  a=@rd
  %+  turn  rds
  |=  b=@rd
  [%rd %n a b (~(add rd %n) a b) (~(sub rd %n) a b) (~(mul rd %n) a b) (~(div rd %n) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rds
  |=  a=@rd
  %+  turn  rds
  |=  b=@rd
  [%rd %u a b (~(add rd %u) a b) (~(sub rd %u) a b) (~(mul rd %u) a b) (~(div rd %u) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rds
  |=  a=@rd
  %+  turn  rds
  |=  b=@rd
  [%rd %d a b (~(add rd %d) a b) (~(sub rd %d) a b) (~(mul rd %d) a b) (~(div rd %d) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rds
  |=  a=@rd
  %+  turn  rds
  |=  b=@rd
  [%rd %z a b (~(add rd %z) a b) (~(sub rd %z) a b) (~(mul rd %z) a b) (~(div rd %z) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rqs
  |=  a=@rq
  %+  turn  rqs
  |=  b=@rq
  [%rq %n a b (~(add rq %n) a b) (~(sub rq %n) a b) (~(mul rq %n) a b) (~(div rq %n) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rqs
  |=  a=@rq
  %+  turn  rqs
  |=  b=@rq
  [%rq %u a b (~(add rq %u) a b) (~(sub rq %u) a b) (~(mul rq %u) a b) (~(div rq %u) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rqs
  |=  a=@rq
  %+  turn  rqs
  |=  b=@rq
  [%rq %d a b (~(add rq %d) a b) (~(sub rq %d) a b) (~(mul rq %d) a b) (~(div rq %d) a b)]
  ^-  (list [@tas @tas @ @ @ @ @ @])
  %-  zing
  %+  turn  rqs
  |=  a=@rq
  %+  turn  rqs
  |=  b=@rq
  [%rq %z a b (~(add rq %z) a b) (~(sub rq %z) a b) (~(mul rq %z) a b) (~(div rq %z) a b)]
==
