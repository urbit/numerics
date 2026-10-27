  ::  /sur/lagoon
::::  Types for Lagoon compatibility
::
|%
+$  ray               ::  $ray:  n-dimensional array
  $:  =meta           ::  descriptor
      data=@ux        ::  data, row-major order, 1-pin MSB
  ==
::
+$  prec  [a=@ b=@]   ::  fixed-point precision, a+b+1=bloq
+$  meta              ::  $meta:  metadata for a $ray
  $:  shape=(list @)  ::  list of dimension lengths
      =bloq           ::  logarithm of bitwidth.  For `%cplx`, bloq is the log₂
                      ::  of the packed-pair element width (both components
                      ::  combined) — e.g., `@cs` (two `@rs` components) has
                      ::  bloq=6, while each component `@rs` has bloq=5.
      =kind           ::  name of data type
      tail=*          ::  per-kind specialization data (~ unless a kind needs
                      ::  it).  Deliberately last and untyped so a kind can add
                      ::  data without rewriting the jet axes.  %fixp uses it for
                      ::  the [a=@ b=@] fixed-point precision; other kinds: ~.
  ==
::
+$  kind              ::  $kind:  type of array scalars
  $?  %i754           ::  IEEE 754 float
      %uint           ::  unsigned integer
      %int2           ::  2s-complement integer (/lib/twoc)
      %unum           ::  unum/posit (/lib/unum) @rpb @rph @rps @rpd
      %cplx           ::  complex (/lib/complex) @ch/@cs/@cd/@cq; bloq=total width
      %fixp           ::  fixed-point Q a.b (/lib/fixed); prec [a b] in meta.tail
  ==
::
+$  baum              ::  $baum:  ndray with metadata
  $:  =meta           ::
      data=ndray      ::
  ==
::
+$  ndray             ::  $ndray:  n-dim array as nested list
    $@  @             ::  single item
    (list ndray)      ::  nonempty list of children, in row-major order
::
+$  slice  (unit [(unit @) (unit @)])
--
