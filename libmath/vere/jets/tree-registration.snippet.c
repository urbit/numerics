// math.hoon transcendental jets -- tree.c registration.
//
// Add this block to EACH per-hoon-kelvin tree in the vere runtime:
//   pkg/noun/jets/135/tree.c, .../136/tree.c, .../137/tree.c
// (the active tree is chosen by the ship's `hoon-version`; numerics ships are
//  136, the stock urbit/urbit ships used for testing are 135).
//
// Place the three static decls just BEFORE `static u3j_core _13X_non_d[]`
// (next to the lagoon `_13X_non__la_core_d` block), then add the "math" entry
// to `_13X_non_d[]` as a sibling of "lagoon".  Replace 13X with 135/136/137.

//  --- decls (before _13X_non_d) ------------------------------------------
static u3j_harm _13X_non__math_rd_exp_a[] = {{".2", u3wi_rd_exp}, {}};
static u3j_core _13X_non__math_rd_d[] =
  { { "exp", 7, _13X_non__math_rd_exp_a, 0, no_hashes },
    {}
  };
static u3j_core _13X_non__math_d[] =
  { { "rd", 7, 0, _13X_non__math_rd_d, no_hashes },
    {}
  };

//  --- entry (inside _13X_non_d[], next to "lagoon") ----------------------
//    { "lagoon", 7, 0, _13X_non__la_core_d, no_hashes },
      { "math",   7, 0, _13X_non__math_d,    no_hashes },
//    { "mice", 7, _13X_non__mice_a, 0, no_hashes },   // (136/135 only)
//    {}
//  };

//  --- header decls (pkg/noun/jets/q.h and w.h) ---------------------------
//  q.h:  u3_noun u3qi_rd_exp(u3_atom);
//  w.h:  u3_noun u3wi_rd_exp(u3_noun);

//  --- build (pkg/noun/build.zig, in the jet .c source list) --------------
//  "jets/i/math.c",
