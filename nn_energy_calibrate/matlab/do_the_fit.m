function do_the_fit( score_terms_to_fit, delG_NN, delG_NN_err, data, seq_labels, coeff_matrix );

prediction_scale_factor = 0.616; % convert to kcal/mol
fprintf(  'Multiplying input scores by %f!\n', prediction_scale_factor );
data.scores = data.scores * 0.616;

if ~exist( 'coeff_matrix' ); coeff_matrix = eye( length( delG_NN ) ); end;
  
fit_idx = get_fit_idx( score_terms_to_fit, data.score_labels );
n_weights = length( fit_idx );

%regularization_stdev = 0.02;
regularization_stdev = 0.2;
regularization_weights = ones( n_weights, 1 ) / ( regularization_stdev^2 );
regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'intermol' ) ) ) = 1/0.01^2;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'fa_rep' ) ) )  =  1/0.01^2;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'rna_sugar_close' ) ) )  =  1/0.01^2;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded' ) ) ) = 0.0;
regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_a' ) ) ) = 0.0;
regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_c' ) ) ) = 0.0;
regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_g' ) ) ) = 0.0;
regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_u' ) ) ) = 0.0;

D = data.scores(:,fit_idx);
data_weights = (1 ./ delG_NN_err.^2);

S = coeff_matrix * D;
A = S' * diag( data_weights ) * S + diag( regularization_weights );
B = S' * diag( data_weights ) * delG_NN + regularization_weights;
new_weights = inv(A) * B;

output_weights( new_weights, data.score_labels( fit_idx ), 0.0 );

delG_fit =  coeff_matrix * data.scores( :, fit_idx ) * new_weights;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots
%%%%%%%%%%%%%%%%%%%%%%%%%%

minE = min( [delG_NN;delG_fit] ) - 1;
maxE = max( [delG_NN;delG_fit] ) + 1;

plot( [minE maxE], [minE maxE], 'k-', 'linew',2 );
hold on

refit_color = [1 0.7 0.7];
plot( delG_fit, delG_NN, '.','color',refit_color );
fprintf( 'fit   rmsd: %8.3f\n', weighted_rmsd( delG_fit, delG_NN, delG_NN_err ) ); 

raw_color = [0.7 0.7 1];

%delG_raw = coeff_matrix * data.scores(:,fit_idx);
delG_raw = coeff_matrix * data.scores(:,1);

plot( delG_raw, delG_NN, '.', 'color', raw_color ); hold on
fprintf( 'raw   rmsd: %8.3f\n', weighted_rmsd( delG_raw, delG_NN, delG_NN_err ) );

% my scratch pad
score = data.scores(:,1);
fa_atr = data.scores( :,find( strcmp( data.score_labels, 'fa_atr' ) ) );
stack_elec = data.scores( :,find( strcmp( data.score_labels, 'stack_elec' ) ) );
hbond_sc = data.scores(:,find( strcmp( data.score_labels, 'hbond_sc' ) ) );
fa_stack = data.scores( :,find( strcmp( data.score_labels, 'fa_stack' ) ) );
%unfolded = data.scores( :,find( strcmp( data.score_labels, 'unfolded' ) ) );
intermol = data.scores( :,find( strcmp( data.score_labels, 'intermol' ) ) );
geom_sol = data.scores( :,find( strcmp( data.score_labels, 'geom_sol' ) ) );
rna_sugar_close = data.scores( :,find( strcmp( data.score_labels, 'rna_sugar_close' ) ) );
unfolded_a = data.scores( :, find( strcmp(   data.score_labels, 'unfolded_a' ) ) );
unfolded_g = data.scores( :, find( strcmp(   data.score_labels, 'unfolded_g' ) ) );
unfolded_u = data.scores( :, find( strcmp(   data.score_labels, 'unfolded_u' ) ) );
unfolded_c = data.scores( :, find( strcmp(   data.score_labels, 'unfolded_c' ) ) );

%delG_refit =  coeff_matrix * (score  + 0.3 * stack_elec );
%plot( delG_refit, delG_NN, '.', 'color', [0.7 1 0.7] ); hold on
%fprintf( 'refit rmsd: %8.3f\n', weighted_rmsd( delG_refit, delG_NN, delG_NN_err ) );

for i = 1:length( delG_fit )
  h = text( delG_raw(i), delG_NN(i), seq_labels{i} );
  set(h,'fontname','courier','interp','none');

  refit_color = fade_color( [1 0 0], 1.5*min( delG_NN_err(i)-0.1, 0.5) );
  raw_color   = fade_color( [0 0 1], 1.5*min( delG_NN_err(i)-0.1, 0.5) );

  plot( delG_fit(i), delG_NN(i), '.', 'color', refit_color );  
  plot( delG_raw(i), delG_NN(i), '.', 'color', raw_color );  

  plot( delG_fit(i) * [1 1], delG_NN(i) + delG_NN_err(i) * [-1 1], '-', 'color', refit_color );  
  plot( delG_raw(i) * [1 1], delG_NN(i) + delG_NN_err(i) * [-1 1], '-','color',raw_color );  
end

hold off
xlabel( 'Rosetta prediction' );
ylabel( 'Turner-rule-based (kcal/mol)' );
legend( '','reweighted','raw',2);
axis( [ minE maxE minE maxE ] );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rmsd =  weighted_rmsd( delG_fit, delG_NN, delG_NN_err )

data_weights = 1./delG_NN_err.^2;
rmsd = sqrt(sum( data_weights' .* (delG_fit - delG_NN ).^2 )/sum(data_weights)) ;