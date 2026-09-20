  ::
::::  /lib/maroon
::
::  MAchine leaRning in hOON: models over Lagoon arrays.
::
::  This is the classical-ML line of urbit/numerics#86 -- models built out of
::  whole-array Lagoon operations and Saloon's decompositions, NOT the tinygrad
::  autodiff port (that is /lib/tinygrad, retargeted separately in #85).
::
::  DATA LAYOUT.  A dataset is an `n x d` %i754 ray: one row per sample, one
::  column per feature.  Labels and indices are %uint bloq-6 rays of shape
::  ~[n], as Lagoon's +argmin-dim returns.  Centroids and components are rays
::  whose rows/columns follow the same convention, so everything composes
::  without reshaping.
::
::  NO PER-ELEMENT LOOPS.  Every arm here is written in terms of whole-array
::  Lagoon ops (+cdist-sq, +argmin-dim, +mean-dim, +mmul, +broadcast-to), so
::  the jets those arms grow will carry these models with them.  The one
::  exception is +update, which must gather rows per cluster.
::
/-  ls=lagoon
/+  *lagoon, saloon, rand, i754rand
::                                                    ::
::::                    ++ma                          ::  machine learning
~%  %maroon  ..part  ~
|%
::    +make:  [rounding-mode @r] -> _ma
::
::  A copy of the +ma core with the rounding mode and tolerance set, as
::  Saloon's +sake does.  Prefer it to bare +ma, whose rtol default is a
::  denormal rather than a usable tolerance.
::  Source
++  make
  |=  [inrnd=rounding-mode inrtol=@r]
  %*(. ma rnd inrnd, rtol inrtol)
::
++  ma
  ^|
  =+  [rnd=`rounding-mode`%n rtol=`@r`0x1]
  ~/  %ma-core
  |%
  ::  the Lagoon and Saloon cores at our settings
  ++  lag  (lake rnd)
  ++  sal  (sake:saloon rnd rtol)
  ::
  ::    +rows / +cols:  the sample and feature counts of an n x d dataset.
  ++  rows  |=(x=ray:ls ^-(@ (snag 0 shape.meta.x)))
  ++  cols  |=(x=ray:ls ^-(@ (snag 1 shape.meta.x)))
  ::    +row:  [x=ray i=@ud] -> ray
  ::
  ::  Row .i of a 2-D ray as a rank-1 ray of shape ~[d].
  ::  Source
  ++  row
    |=  [x=ray:ls i=@ud]
    ^-  ray:ls
    =/  d  (cols x)
    =/  o  (zeros:lag [~[d] bloq.meta.x %i754 ~])
    =/  j  0
    |-  ^-  ray:ls
    ?:  =(j d)  o
    $(j +(j), o (set-item:lag o ~[j] (get-item:lag x ~[i j])))
  ::    +set-row:  [x=ray i=@ud v=ray] -> ray
  ::
  ::  .x with row .i replaced by the rank-1 ray .v.
  ::  Source
  ++  set-row
    |=  [x=ray:ls i=@ud v=ray:ls]
    ^-  ray:ls
    =/  d  (cols x)
    =/  j  0
    |-  ^-  ray:ls
    ?:  =(j d)  x
    $(j +(j), x (set-item:lag x ~[i j] (get-item:lag v ~[j])))
  ::    +col-mean:  x=ray -> ray
  ::
  ::  The mean of each feature (column) of an n x d dataset, as a rank-1 ray
  ::  of shape ~[d].
  ::  Source
  ++  col-mean  |=(x=ray:ls ^-(ray:ls (mean-dim:lag x 0)))
  ::    +center:  x=ray -> ray
  ::
  ::  .x with each column's mean subtracted, the standard first step of PCA.
  ::  Source
  ++  center
    |=  x=ray:ls
    ^-  ray:ls
    (sub:lag x (broadcast-to:lag (col-mean x) shape.meta.x))
  ::    +cov:  [x=ray ddof=@ud] -> ray
  ::
  ::  The d x d covariance matrix of .x: Xc^T*Xc / (n - ddof), with .ddof=1 the
  ::  usual sample covariance.  +pca does not use this (it decomposes Xc
  ::  directly, which is better conditioned), but it is what a caller wants for
  ::  a Gaussian model or a whitening transform.
  ::  Source
  ++  cov
    |=  [x=ray:ls ddof=@ud]
    ^-  ray:ls
    =/  n  (rows x)
    ~|  'maroon cov: ddof must be less than the sample count'
    ?>  (gth n ddof)
    =/  xc  (center x)
    %+  div-scalar:lag
      (mmul:lag (transpose:lag xc) xc)
    (i754-sun:lag bloq.meta.x (sub n ddof))
  ::    +pca:  [x=ray k=@ud] -> [comp=ray vals=ray mean=ray]
  ::
  ::  Principal component analysis of an n x d dataset, keeping .k components:
  ::  .comp is d x k (the components as COLUMNS, largest variance first),
  ::  .vals is a rank-1 ray of the variance each explains, and .mean is the
  ::  feature mean +pca-transform needs.
  ::
  ::  Computed from the SVD of the centered data rather than an eigen-
  ::  decomposition of the covariance: forming Xc^T*Xc squares the condition
  ::  number, and Saloon's one-sided Jacobi already returns its singular values
  ::  in descending order, which is the order PCA wants.  The explained
  ::  variances are s_i^2/(n-1).
  ::
  ::  Component SIGNS are arbitrary (as in every PCA implementation): a
  ::  component and its negation describe the same axis.  Compare with a sign
  ::  convention or on absolute values.
  ::    Examples
  ::      > =ma  (make %n .~1e-12)
  ::      > =x  (en-ray:la [[~[3 2] 6 %i754 ~] ~[~[.~1 .~0] ~[.~2 .~0] ~[.~3 .~0]]])
  ::      > ;;((list @rd) data:(de-ray:la comp:(pca:ma x 1)))
  ::      ~[.~1 .~0]                       ::  the x-axis, exactly
  ::      > ;;((list @rd) data:(de-ray:la vals:(pca:ma x 1)))
  ::      ~[.~0.9999999999999998]          ::  variance 1, via sqrt-then-square
  ::  Source
  ++  pca
    |=  [x=ray:ls k=@ud]
    ^-  [comp=ray:ls vals=ray:ls mean=ray:ls]
    ?>  =(2 (lent shape.meta.x))
    =/  n  (rows x)
    =/  d  (cols x)
    ~|  'maroon pca: k must be in 1..d'
    ?>  ?&(!=(0 k) (lte k d))
    ~|  'maroon pca: needs at least two samples'
    ?>  (gth n 1)
    =/  xc   (center x)
    =/  res  (svd:sal xc)
    ::  the top k right-singular vectors, as columns
    =/  comp  (submatrix:lag ~[~ `[`0 `(dec k)]] v.res)
    ::  variance explained: s_i^2 / (n-1)
    =/  sk  (submatrix:lag ~[`[`0 `(dec k)]] s.res)
    =/  vals
      %+  div-scalar:lag  (mul:lag sk sk)
      (i754-sun:lag bloq.meta.x (dec n))
    [comp vals (col-mean x)]
  ::    +pca-transform:  [x=ray comp=ray mean=ray] -> ray
  ::
  ::  Projects .x (n x d) onto the components: subtract .mean, then multiply by
  ::  .comp, giving n x k scores.  Takes .comp/.mean rather than re-fitting, so
  ::  new data lands in the same basis as the training set.
  ::  Source
  ++  pca-transform
    |=  [x=ray:ls comp=ray:ls mean=ray:ls]
    ^-  ray:ls
    (mmul:lag (sub:lag x (broadcast-to:lag mean shape.meta.x)) comp)
  ::
  ::  k-means.  Lloyd's algorithm is deterministic given its starting
  ::  centroids, so it is kept separate from the (random) k-means++ seeding:
  ::  +kmeans takes the initial centroids, +kmeans-pp produces them, and
  ::  +kmeans-fit is the two together.
  ::
  ::    +assign:  [x=ray c=ray] -> ray
  ::
  ::  The index of the nearest centroid for each sample, as a %uint bloq-6 ray
  ::  of shape ~[n].  One +cdist-sq and one +argmin-dim, so the whole
  ::  assignment step is two Lagoon calls; ties go to the lower centroid index.
  ::  Source
  ++  assign
    |=  [x=ray:ls c=ray:ls]
    ^-  ray:ls
    (argmin-dim:lag (cdist-sq:lag x c) 1)
  ::    +update:  [x=ray lab=ray c=ray k=@ud] -> ray
  ::
  ::  New centroids: the mean of the samples assigned to each cluster.  An
  ::  EMPTY cluster keeps its previous centroid (the common convention, and it
  ::  keeps the shape fixed) rather than being dropped or reseeded.
  ::  Source
  ++  update
    |=  [x=ray:ls lab=ray:ls c=ray:ls k=@ud]
    ^-  ray:ls
    =/  n  (rows x)
    =/  d  (cols x)
    =/  b  bloq.meta.x
    =/  out  c
    =/  j  0
    |-  ^-  ray:ls
    ?:  =(j k)  out
    ::  sum the rows of cluster j, counting them
    =/  acc
      =/  i    0
      =/  cnt  0
      =/  sum  (zeros:lag [~[d] b %i754 ~])
      |-  ^-  [cnt=@ud sum=ray:ls]
      ?:  =(i n)  [cnt sum]
      ?.  =(j `@`(get-item:lag lab ~[i]))
        $(i +(i))
      $(i +(i), cnt +(cnt), sum (add:lag sum (row x i)))
    ?:  =(0 cnt.acc)  $(j +(j))
    %=  $
      j    +(j)
      out  (set-row out j (div-scalar:lag sum.acc (i754-sun:lag b cnt.acc)))
    ==
  ::    +kmeans:  [x=ray c0=ray maxit=@ud] -> [labels=ray centroids=ray iter=@ud]
  ::
  ::  Lloyd's algorithm from the given initial centroids .c0 (k x d).  Stops
  ::  when the assignment stops changing -- the exact fixed point, not a
  ::  tolerance -- or after .maxit iterations, reporting the count so a caller
  ::  can tell convergence from exhaustion.  Deterministic: same input, same
  ::  output, every time.
  ::  Source
  ++  kmeans
    |=  [x=ray:ls c0=ray:ls maxit=@ud]
    ^-  [labels=ray:ls centroids=ray:ls iter=@ud]
    ?>  =(2 (lent shape.meta.x))
    ?>  =(2 (lent shape.meta.c0))
    ~|  'maroon kmeans: centroids and samples must agree on the feature count'
    ?>  =((cols x) (cols c0))
    =/  k    (rows c0)
    =/  c    c0
    =/  lab  (assign x c)
    =/  it   0
    |-  ^-  [labels=ray:ls centroids=ray:ls iter=@ud]
    ?:  =(it maxit)  [lab c it]
    =/  c2    (update x lab c k)
    =/  lab2  (assign x c2)
    ?:  =(lab lab2)  [lab2 c2 +(it)]
    $(it +(it), c c2, lab lab2)
  ::    +kmeans-pp:  [x=ray k=@ud r=rng] -> [centroids=ray r=rng]
  ::
  ::  k-means++ seeding (Arthur & Vassilvitskii 2007): the first centroid is a
  ::  uniformly drawn sample, and each subsequent one is drawn with probability
  ::  proportional to its squared distance from the centroids chosen so far, so
  ::  the seeds spread out instead of clumping.  Threads the RNG through
  ::  explicitly, so a given seed always gives the same centroids.
  ::
  ::  The weighted draw goes through /lib/i754rand's alias table, which is
  ::  built per step (O(n)) and drawn from once.
  ::  Source
  ++  kmeans-pp
    |=  [x=ray:ls k=@ud r=rng:rand]
    ^-  [centroids=ray:ls r=rng:rand]
    ?>  =(2 (lent shape.meta.x))
    =/  n  (rows x)
    =/  d  (cols x)
    =/  b  bloq.meta.x
    ~|  'maroon kmeans-pp: k must be in 1..n'
    ?>  ?&(!=(0 k) (lte k n))
    ::  first centroid: a uniformly chosen sample
    =^  i0  r  (below:uni:rand r n)
    =/  c  (set-row (zeros:lag [~[k d] b %i754 ~]) 0 (row x i0))
    =/  j  1
    |-  ^-  [centroids=ray:ls r=rng:rand]
    ?:  =(j k)  [c r]
    ::  squared distance from each sample to its nearest chosen centroid
    =/  chosen  (submatrix:lag ~[`[`0 `(dec j)] ~] c)
    =/  dist2   (min-dim:lag (cdist-sq:lag x chosen) 1)
    ::  draw index j with probability proportional to that distance
    =/  weights=(list @rd)
      %+  turn  (gulf 0 (dec n))
      |=(i=@ `@rd`(change-bloq b (get-item:lag dist2 ~[i])))
    =^  pick  r  (categorical:rd:dist:i754rand r (build:alias:i754rand weights))
    $(j +(j), c (set-row c j (row x pick)))
  ::    +change-bloq:  [b=@ v=@] -> @rd
  ::
  ::  One %i754 scalar of width .b as an @rd, which is the width
  ::  /lib/i754rand's alias table fixes its weights at.
  ::  Source
  ++  change-bloq
    |=  [b=@ v=@]
    ^-  @rd
    ?:  =(6 b)  `@rd`v
    =/  one  (scalar-to-ray:lag [~[1] b %i754 ~] v)
    `@rd`(get-item:lag (change:lag one %i754 6) ~[0])
  ::    +kmeans-fit:  [x=ray k=@ud maxit=@ud r=rng]
  ::                  -> [labels=ray centroids=ray iter=@ud r=rng]
  ::
  ::  k-means++ seeding followed by Lloyd's algorithm.
  ::    Examples
  ::      > =ma  (make %n .~1e-12)
  ::      > =r   (from-atom:seed:rand %sm64 0)
  ::      > =x   (en-ray:la [[~[4 1] 6 %i754 ~] ~[.~0 .~1 .~10 .~11]])
  ::      > (sort ;;((list @) data:(de-ray:la labels:(kmeans-fit:ma x 2 20 r))) lth)
  ::      ~[0 0 1 1]
  ::  Source
  ++  kmeans-fit
    |=  [x=ray:ls k=@ud maxit=@ud r=rng:rand]
    ^-  [labels=ray:ls centroids=ray:ls iter=@ud r=rng:rand]
    =/  seed  (kmeans-pp x k r)
    =/  res   (kmeans x centroids.seed maxit)
    [labels.res centroids.res iter.res r.seed]
  ::    +inertia:  [x=ray lab=ray c=ray] -> @
  ::
  ::  The k-means objective: the sum of squared distances from each sample to
  ::  its assigned centroid, summed left to right.  Lloyd's algorithm decreases
  ::  it monotonically, so it is the number to compare two runs with.
  ::  Source
  ++  inertia
    |=  [x=ray:ls lab=ray:ls c=ray:ls]
    ^-  @
    =/  n  (rows x)
    =/  b  bloq.meta.x
    =/  d2  (cdist-sq:lag x c)
    =/  i    0
    =/  acc  (i754-sun:lag b 0)
    |-  ^-  @
    ?:  =(i n)  acc
    =/  j  `@`(get-item:lag lab ~[i])
    $(i +(i), acc (fadd:sal b acc (get-item:lag d2 ~[i j])))
  --
--
