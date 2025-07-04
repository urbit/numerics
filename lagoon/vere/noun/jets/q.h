/// @file

#ifndef U3_JETS_Q_H
#define U3_JETS_Q_H

#include "types.h"

  /** Tier 1.
  **/
    u3_noun u3qa_add(u3_atom, u3_atom);
    u3_noun u3qa_dec(u3_atom);
    u3_noun u3qa_div(u3_atom, u3_atom);
    u3_noun u3qa_gte(u3_atom, u3_atom);
    u3_noun u3qa_gth(u3_atom, u3_atom);
    u3_noun u3qa_inc(u3_atom);
    u3_noun u3qa_lte(u3_atom, u3_atom);
    u3_noun u3qa_lth(u3_atom, u3_atom);
    u3_noun u3qa_max(u3_atom, u3_atom);
    u3_noun u3qa_min(u3_atom, u3_atom);
    u3_noun u3qa_mod(u3_atom, u3_atom);
    u3_noun u3qa_mul(u3_atom, u3_atom);
    u3_noun u3qa_sub(u3_atom, u3_atom);

  /** Tier 2.
  **/
    u3_noun u3qb_bind(u3_noun, u3_noun);
    u3_noun u3qb_clap(u3_noun, u3_noun, u3_noun);
    u3_noun u3qb_drop(u3_noun);
    u3_noun u3qb_flop(u3_noun);
    u3_noun u3qb_lent(u3_noun);
    u3_noun u3qb_levy(u3_noun, u3_noun);
    u3_noun u3qb_lien(u3_noun, u3_noun);
    u3_noun u3qb_murn(u3_noun, u3_noun);
    u3_noun u3qb_need(u3_noun);
    u3_noun u3qb_reap(u3_atom, u3_noun);
    u3_noun u3qb_reel(u3_noun, u3_noun);
    u3_noun u3qb_roll(u3_noun, u3_noun);
    u3_noun u3qb_skid(u3_noun, u3_noun);
    u3_noun u3qb_skim(u3_noun, u3_noun);
    u3_noun u3qb_skip(u3_noun, u3_noun);
    u3_noun u3qb_scag(u3_atom, u3_noun);
    u3_noun u3qb_slag(u3_atom, u3_noun);
    u3_noun u3qb_snag(u3_atom, u3_noun);
    u3_noun u3qb_sort(u3_noun, u3_noun);
    u3_noun u3qb_turn(u3_noun, u3_noun);
    u3_noun u3qb_weld(u3_noun, u3_noun);

  /** Tier 3.
  **/
    u3_noun u3qc_bex(u3_atom);
    u3_noun u3qc_xeb(u3_atom);
    u3_noun u3qc_can(u3_atom, u3_noun);
    u3_noun u3qc_cap(u3_atom);
    u3_noun u3qc_cat(u3_atom, u3_atom, u3_atom);
    u3_noun u3qc_con(u3_atom, u3_atom);
    u3_noun u3qc_cut(u3_atom, u3_atom, u3_atom, u3_atom);
    u3_noun u3qc_dis(u3_atom, u3_atom);
    u3_noun u3qc_dor(u3_atom, u3_atom);
    u3_noun u3qc_dvr(u3_atom, u3_atom);
    u3_noun u3qc_end(u3_atom, u3_atom, u3_atom);
    u3_noun u3qc_gor(u3_atom, u3_atom);
    u3_noun u3qc_lsh(u3_atom, u3_atom, u3_atom);
    u3_noun u3qc_mas(u3_atom);
    u3_noun u3qc_met(u3_atom, u3_atom);
    u3_noun u3qc_mix(u3_atom, u3_atom);
    u3_noun u3qc_mor(u3_atom, u3_atom);
    u3_noun u3qc_muk(u3_atom, u3_atom, u3_atom);
    u3_noun u3qc_peg(u3_atom, u3_atom);
    u3_noun u3qc_pow(u3_atom, u3_atom);
    u3_noun u3qc_rap(u3_atom, u3_noun);
    u3_noun u3qc_rep(u3_atom, u3_atom, u3_noun);
    u3_noun u3qc_rev(u3_atom, u3_atom, u3_atom);
    u3_noun u3qc_rip(u3_atom, u3_atom, u3_atom);
    u3_noun u3qc_rsh(u3_atom, u3_atom, u3_atom);
    u3_noun u3qc_swp(u3_atom, u3_atom);
    u3_noun u3qc_sqt(u3_atom);

    u3_noun u3_po_find_prefix(c3_y one, c3_y two, c3_y three);
    u3_noun u3_po_find_suffix(c3_y one, c3_y two, c3_y three);
    void u3_po_to_prefix(u3_noun id, c3_y* a, c3_y* b, c3_y* c);
    void u3_po_to_suffix(u3_noun id, c3_y* a, c3_y* b, c3_y* c);

  /** Tier 4.
  **/
    u3_noun u3qdb_all(u3_noun, u3_noun);
    u3_noun u3qdb_any(u3_noun, u3_noun);
    u3_noun u3qdb_apt(u3_noun);
    u3_noun u3qdb_bif(u3_noun, u3_noun);
    u3_noun u3qdb_dif(u3_noun, u3_noun);
    u3_noun u3qdb_gas(u3_noun, u3_noun);
    u3_noun u3qdb_get(u3_noun, u3_noun);
    u3_noun u3qdb_has(u3_noun, u3_noun);
    u3_noun u3qdb_int(u3_noun, u3_noun);
    u3_noun u3qdb_key(u3_noun);
    u3_noun u3qdb_put(u3_noun, u3_noun, u3_noun);
    u3_noun u3qdb_run(u3_noun, u3_noun);
#   define u3qdb_tap u3qdi_tap
    u3_noun u3qdb_uni(u3_noun, u3_noun);
    u3_noun u3qdb_urn(u3_noun, u3_noun);
#   define u3qdb_wyt u3qdi_wyt

    u3_noun u3qdi_apt(u3_noun);
    u3_noun u3qdi_bif(u3_noun, u3_noun);
    u3_noun u3qdi_dif(u3_noun, u3_noun);
    u3_noun u3qdi_gas(u3_noun, u3_noun);
    u3_noun u3qdi_has(u3_noun, u3_noun);
    u3_noun u3qdi_int(u3_noun, u3_noun);
    u3_noun u3qdi_put(u3_noun, u3_noun);
    u3_noun u3qdi_rep(u3_noun, u3_noun);
    u3_noun u3qdi_run(u3_noun, u3_noun);
    u3_noun u3qdi_tap(u3_noun);
    u3_noun u3qdi_uni(u3_noun, u3_noun);
    u3_noun u3qdi_wyt(u3_noun);

  /** Tier 5.
  **/
    u3_noun u3qe_cue(u3_atom);
    u3_noun u3qe_jam(u3_atom);
    u3_noun u3qe_mat(u3_atom);
    u3_noun u3qe_rub(u3_atom, u3_atom);
    u3_noun u3qe_leer(u3_atom);
    u3_noun u3qe_lore(u3_atom);
    u3_noun u3qe_loss(u3_noun, u3_noun);
    u3_noun u3qe_lune(u3_atom);
    u3_noun u3qe_repg(u3_noun, u3_noun, u3_noun);
    u3_noun u3qe_rexp(u3_noun, u3_noun);
    u3_noun u3qe_trip(u3_atom);

    u3_atom u3qe_scot(u3_atom, u3_atom);
    u3_atom u3qe_scow(u3_atom, u3_atom);

    u3_noun u3qea_ecba_en(u3_atom, u3_atom);
    u3_noun u3qea_ecba_de(u3_atom, u3_atom);
    u3_noun u3qea_ecbb_en(u3_atom, u3_atom);
    u3_noun u3qea_ecbb_de(u3_atom, u3_atom);
    u3_noun u3qea_ecbc_en(u3_atom, u3_atom);
    u3_noun u3qea_ecbc_de(u3_atom, u3_atom);

    u3_noun u3qea_cbca_en(u3_atom, u3_atom, u3_atom);
    u3_noun u3qea_cbca_de(u3_atom, u3_atom, u3_atom);
    u3_noun u3qea_cbcb_en(u3_atom, u3_atom, u3_atom);
    u3_noun u3qea_cbcb_de(u3_atom, u3_atom, u3_atom);
    u3_noun u3qea_cbcc_en(u3_atom, u3_atom, u3_atom);
    u3_noun u3qea_cbcc_de(u3_atom, u3_atom, u3_atom);

    u3_noun u3qea_de(u3_atom, u3_atom);
    u3_noun u3qea_en(u3_atom, u3_atom);

    u3_atom u3qe_fein_ob(u3_atom pyn);
    u3_atom u3qe_fynd_ob(u3_atom pyn);

    u3_noun u3qe_hmac(u3_noun, u3_atom, u3_atom,
                      u3_atom, u3_atom, u3_atom, u3_atom);

    u3_noun u3qe_en_base16(u3_atom len, u3_atom dat);
    u3_noun u3qe_de_base16(u3_atom inp);

    u3_noun u3qe_json_de(u3_atom);
    u3_atom u3qe_json_en(u3_noun);

    u3_noun u3qeo_raw(u3_atom, u3_atom);

    u3_noun u3qef_drg(u3_noun, u3_atom);
    u3_noun u3qef_lug(u3_noun, u3_noun, u3_atom, u3_atom);

    u3_noun u3qer_add(u3_atom, u3_atom, u3_atom);
    u3_noun u3qer_sub(u3_atom, u3_atom, u3_atom);
    u3_noun u3qer_mul(u3_atom, u3_atom, u3_atom);
    u3_noun u3qer_div(u3_atom, u3_atom, u3_atom);
    u3_noun u3qer_sqt(u3_atom, u3_atom);
    u3_noun u3qer_fma(u3_atom, u3_atom, u3_atom, u3_atom);
    u3_noun u3qer_lth(u3_atom, u3_atom);
    u3_noun u3qer_lte(u3_atom, u3_atom);
    u3_noun u3qer_equ(u3_atom, u3_atom);
    u3_noun u3qer_gte(u3_atom, u3_atom);
    u3_noun u3qer_gth(u3_atom, u3_atom);

    u3_noun u3qet_add(u3_atom, u3_atom, u3_atom);
    u3_noun u3qet_sub(u3_atom, u3_atom, u3_atom);
    u3_noun u3qet_mul(u3_atom, u3_atom, u3_atom);
    u3_noun u3qet_div(u3_atom, u3_atom, u3_atom);
    u3_noun u3qet_sqt(u3_atom, u3_atom);
    u3_noun u3qet_fma(u3_atom, u3_atom, u3_atom, u3_atom);
    u3_noun u3qet_lth(u3_atom, u3_atom);
    u3_noun u3qet_lte(u3_atom, u3_atom);
    u3_noun u3qet_equ(u3_atom, u3_atom);
    u3_noun u3qet_gte(u3_atom, u3_atom);
    u3_noun u3qet_gth(u3_atom, u3_atom);

    u3_noun u3qeq_add(u3_atom, u3_atom, u3_atom);
    u3_noun u3qeq_sub(u3_atom, u3_atom, u3_atom);
    u3_noun u3qeq_mul(u3_atom, u3_atom, u3_atom);
    u3_noun u3qeq_div(u3_atom, u3_atom, u3_atom);
    u3_noun u3qeq_sqt(u3_atom, u3_atom);
    u3_noun u3qeq_fma(u3_atom, u3_atom, u3_atom, u3_atom);
    u3_noun u3qeq_lth(u3_atom, u3_atom);
    u3_noun u3qeq_lte(u3_atom, u3_atom);
    u3_noun u3qeq_equ(u3_atom, u3_atom);
    u3_noun u3qeq_gte(u3_atom, u3_atom);
    u3_noun u3qeq_gth(u3_atom, u3_atom);

    u3_noun u3qes_add(u3_atom, u3_atom, u3_atom);
    u3_noun u3qes_sub(u3_atom, u3_atom, u3_atom);
    u3_noun u3qes_mul(u3_atom, u3_atom, u3_atom);
    u3_noun u3qes_div(u3_atom, u3_atom, u3_atom);
    u3_noun u3qes_sqt(u3_atom, u3_atom);
    u3_noun u3qes_fma(u3_atom, u3_atom, u3_atom, u3_atom);
    u3_noun u3qes_lth(u3_atom, u3_atom);
    u3_noun u3qes_lte(u3_atom, u3_atom);
    u3_noun u3qes_equ(u3_atom, u3_atom);
    u3_noun u3qes_gte(u3_atom, u3_atom);
    u3_noun u3qes_gth(u3_atom, u3_atom);

  /** Tier 6.
  **/
    u3_noun u3qf_bull(u3_noun, u3_noun);
    u3_noun u3qf_cell(u3_noun, u3_noun);
    u3_noun u3qf_comb(u3_noun, u3_noun);
    u3_noun u3qf_cons(u3_noun, u3_noun);
    u3_noun u3qf_core(u3_noun, u3_noun);
    u3_noun u3qf_cube(u3_noun, u3_noun);
    u3_noun u3qf_face(u3_noun, u3_noun);
    u3_noun u3qf_fine(u3_noun, u3_noun, u3_noun);
    u3_noun u3qf_fitz(u3_noun, u3_noun);
    u3_noun u3qf_flan(u3_noun, u3_noun);
    u3_noun u3qf_flay(u3_noun);
    u3_noun u3qf_flip(u3_noun);
    u3_noun u3qf_flor(u3_noun, u3_noun);
    u3_noun u3qf_forq(u3_noun, u3_noun);
    u3_noun u3qf_fork(u3_noun);
    u3_noun u3qf_grof(u3_noun);
    u3_noun u3qf_hint(u3_noun, u3_noun);
    u3_noun u3qf_hike(u3_noun, u3_noun);
    u3_noun u3qf_look(u3_noun, u3_noun);
    u3_noun u3qf_loot(u3_noun, u3_noun);
    u3_noun u3qf_slot(u3_atom, u3_noun);
    u3_noun u3qf_type(u3_noun);

    u3_noun u3qfl_bunt(u3_noun, u3_noun);
    u3_noun u3qfl_whip(u3_noun, u3_noun, u3_noun);

    u3_noun u3qfr_fish(u3_noun, u3_noun, u3_noun, u3_noun);

    u3_noun u3qfp_hack(u3_noun, u3_noun);
    u3_noun u3qfp_late(u3_noun);
    u3_noun u3qfp_open(u3_noun, u3_noun, u3_noun);
    u3_noun u3qfp_nepo(u3_noun, u3_noun);
    u3_noun u3qfp_rake(u3_noun);

    u3_noun u3qi_la_add_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_sub_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_mul_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_div_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_mod_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_adds_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_subs_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_muls_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_divs_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_mods_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_dot_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_diag(u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_transpose(u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_cumsum_i754(u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_argmin_i754(u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_argmax_i754(u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_ravel_i754(u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_min_i754(u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_max_i754(u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_linspace_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_range_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_abs_i754(u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_gth_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_gte_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_lth_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_lte_i754(u3_noun, u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_trace_i754(u3_noun, u3_noun, u3_noun);
    u3_noun u3qi_la_mmul_i754(u3_noun, u3_noun, u3_noun, u3_noun, u3_noun);

#   define u3qfu_van_fan  28
#   define u3qfu_van_rib  58
#   define u3qfu_van_vet  59

    void u3qf_test(const c3_c*, u3_noun);

#endif /* ifndef U3_JETS_Q_H */
