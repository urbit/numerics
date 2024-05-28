  ::  /lib/bmp
::::  Version ~2024.5.24 by ~lagrev-nocfep
::
::    A bitmap file as created by this library consists of
::    24-bit RGB pixel data with a 54-byte header.
::
/-  bmp,
    ls=lagoon
/+  *lagoon
|%
::  LSB-first BMP file header, including width in bytes.
++  fil-hdr
  |=  =meta:ls
  ^-  [dat=@ux len=@]
  =|  val=@ux
  ::  Header field 'BM' (2 bytes)
  =.  val  0x4d42
  ::  File size in bytes (4 bytes)
  =/  syz  (roll shape.meta add)
  =.  val  (con (lsh [3 2] val) syz)
  ::  Reserved (4 bytes)
  =/  res  0x0
  =.  val  (con (lsh [3 4] val) res)
  ::  Start of pixel data (4 bytes) = 54
  =/  off  0x36
  =.  val  (con (lsh [3 4] val) off)
  [val 14]
::  DIB header, including width in bytes.
++  dib-hdr
  |=  =meta:ls
  ^-  [dat=@ux len=@]
  =|  val=@ux
  ::  Header name
  =.  val  `@ux`%'BITMAPINFOHEADER'
  ::  Header size (4 bytes)
  =/  syz  0x28
  =.  val  (con (lsh [3 4] val) syz)
  ::  Image width in pixels (4 bytes)
  =/  wyd  (snag 1 shape.meta)
  =.  val  (con (lsh [3 4] val) wyd)
  ::  Image height in pixels (4 bytes)
  =/  hyt  (snag 0 shape.meta)
  =.  val  (con (lsh [3 4] val) hyt)
  ::  Number of color planes (2 bytes)
  =/  pln  0x1
  =.  val  (con (lsh [3 2] val) pln)
  ::  Bits per pixel (2 bytes)
  =/  bpp  0x18
  =.  val  (con (lsh [3 2] val) bpp)
  ::  Compression method (4 bytes)
  =/  cmp  0x0
  =.  val  (con (lsh [3 4] val) cmp)
  ::  Image size in bytes (4 bytes)
  =/  siz  (roll shape.meta add)
  =.  val  (con (lsh [3 4] val) siz)
  ::  Horizontal resolution in pixels per meter (4 bytes)
  =/  xpm  0x24
  =.  val  (con (lsh [3 4] val) xpm)
  ::  Vertical resolution in pixels per meter (4 bytes)
  =/  ypm  0x24
  =.  val  (con (lsh [3 4] val) ypm)
  ::  Number of colors in the color palette (4 bytes)
  =/  clr  0x0
  =.  val  (con (lsh [3 4] val) clr)
  ::  Number of important colors used (4 bytes)
  =/  imp  0x0
  =.  val  (con (lsh [3 4] val) imp)
  [val 40]
::  ++write-task
::
::  Write to Clay.
::  (write-task /path/to/output/png my-ray)
++  write-task
  |=  [=path =raw:bmp]
  ^-  task:clay
  :*  %info
      %base
      ^-  nori:clay
      :*  %&
      ^-  soba:clay
      :~  :-  path
          ^-  miso:clay
          :-  %ins
          :-  %bmp
          !>((export raw))
  ==  ==  ==
::  ++import
::
::  Import ray from BMP file as RGB.
::  (import /path/to/input/bmp)
++  import
  |=  =path
  ^-  raw:bmp
  =/  fil  .^(@ux %cx path)
  ::  width at 19--22
  =/  wyd  (cut 3 [18 4] fil)
  ::  height at 23--26
  =/  hyt  (cut 3 [22 4] fil)
  ~&  >  [wyd hyt]
  ::  color depth
  =/  dep  (cut 3 [28 4] fil)
  =/  bpm  (div dep 8)
  ::  pixel data at offset
  =/  off  (cut 3 [10 4] fil)
  =/  len  (met 3 fil)
  =/  raw  (cut 3 [off len] fil)
  ~|  %bmp-wrong-shape
  ?>  =(:(mul bpm hyt wyd) (cut 3 [34 4] fil))
  ^-  raw:bmp
  :-  [~[hyt wyd bpm] 3 %uint ~]
  raw
::  ++export
::
::  Export ray to BMP-style hexadecimal.
::  (export my-ray)
++  export
  |=  =raw:bmp
  ^-  bmp:bmp
  |^
  =|  val=@ux
  ::  BMP header
  =.  val  `@ux`-:(fil-hdr meta.raw)
  ::  DIB header
  =.  val  (con (lsh [3 14] val) `@ux`-:(dib-hdr meta.raw))
  ::  Pixel data
  =.  val  `@ux`(con (lsh [3 40] val) data:(unspac:la raw))
  val
  ++  zip3  :: derived but not identical to /lib/sequent
    |*  [a=(list) b=(list) c=(list)]
    =/  d=(list _?>(?=(^ a) i.a))  ~
    |-  ^-  (list _?>(?=(^ a) i.a))
    ?~  a  ?~  b  ?~  c  d
           ~|('lists of unequal length' !!)
           ~|('lists of unequal length' !!)
    ?~  b  ~|('lists of unequal length' !!)
    ?~  c  ~|('lists of unequal length' !!)
    $(a t.a, b t.b, c t.c, d [i.a i.b i.c d])
  --
--
