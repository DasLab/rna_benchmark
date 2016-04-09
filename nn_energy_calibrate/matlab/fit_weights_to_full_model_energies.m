[data,tags] = get_scores_and_tags( '../helix/pdb_scores.out' );
data = add_separate_unfolded_energies( data, tags, '../helix/rna_helix_stack_elec.wts'  ); 

turner_rules = init_delG_NN();
[ seq_labels, delG_NN ] = get_delG_NN_for_tags( tags, turner_rules );
delG_raw = data.scores(:,1);
if isempty( find( strcmp( data.score_labels, 'intermol' ) ) )
  data.score_labels{ end+1 } = 'intermol';
  data.scores( :, end+1 ) = 1.0;
end

% refit weights.
score_terms_to_fit = { 'fa_atr', 'fa_rep',  'fa_stack', 'hbond_sc', 'rna_torsion', 'geom_sol',  'stack_elec' };
%score_terms_to_fit = { 'fa_atr', 'fa_rep',  'fa_stack', 'hbond_sc', 'rna_torsion', 'CI_geom_sol',  'stack_elec', 'intermol', 'unfolded' };
%score_terms_to_fit = { 'fa_atr', 'fa_rep',  'fa_stack', 'hbond_sc', 'rna_torsion', 'geom_sol',  'stack_elec', 'unfolded', 'intermol', 'lk_nonpolar' };
%score_terms_to_fit = { 'fa_atr', 'fa_rep',  'fa_stack', 'hbond_sc', 'rna_torsion', 'occ_sol_fitted',  'stack_elec', 'intermol', 'unfolded' };
%score_terms_to_fit = { 'fa_atr', 'fa_rep',  'fa_stack', 'hbond_sc', 'rna_torsion', 'lk_ball',  'stack_elec', 'intermol', 'unfolded' };
%score_terms_to_fit = { 'score','stack_elec' };
%score_terms_to_fit = { 'score','unfolded_g','unfolded_a','unfolded_c','unfolded_u' };

%score_terms_to_fit = setdiff( data.score_labels, 'score' );
score_terms_to_fit = setdiff( data.score_labels, {'score','ref','unfolded'} );
%score_terms_to_fit = setdiff( data.score_labels, {'score','ref','stack_elec_base_base','stack_elec_base_bb'} );

delG_NN_err = delG_NN*0 + 0.3;
do_the_fit( score_terms_to_fit, delG_NN, delG_NN_err', data, tags );
