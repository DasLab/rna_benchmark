function fit_idx = get_fit_idx( score_terms_to_fit, score_labels );

fit_idx = [];
for i = 1:length(score_labels)
  if ~isempty( find( strcmp( score_terms_to_fit, score_labels{i} ) ) )
    fit_idx = [fit_idx, i ];
  end
end


assert( length( fit_idx ) == length( score_terms_to_fit ) );
