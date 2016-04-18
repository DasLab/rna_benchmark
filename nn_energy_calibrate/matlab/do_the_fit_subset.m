function do_the_fit_subset( score_terms_to_fit, delG_NN, delG_NN_err, data, seq_labels, coeff_matrix );

prediction_scale_factor = 0.4; % convert to kcal/mol
%prediction_scale_factor = 0.616; % convert to kcal/mol if you assume rosetta scores are in kT
fprintf(  'Multiplying input scores by %f!\n', prediction_scale_factor );
data.scores = data.scores * prediction_scale_factor;

% This happens in fit_weights_to_full_model_energies
if ~exist( 'coeff_matrix' ); coeff_matrix = eye( length( delG_NN ) ); end;


%%%%%%%% Fit weights on a subset of score terms and a single weight on the sum of the rest of the scores %%%%%%
%fit_idx = get_fit_idx( score_terms_to_fit, data.score_labels );
%n_weights = length( fit_idx )+1;
%
%regularization_stdev = 0.2;
%regularization_weights = ones( n_weights, 1 ) / ( regularization_stdev^2 );
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'intermol' ) ) ) = 1/0.01^2;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_a' ) ) ) = 0.0;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_c' ) ) ) = 0.0;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_g' ) ) ) = 0.0;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_u' ) ) ) = 0.0;
%
%D = horzcat( data.scores(:,fit_idx), data.scores(:,1) - sum(data.scores(:,fit_idx),2));
%data_weights = (1 ./ delG_NN_err.^2);
%
%S = coeff_matrix * D;
%% with regularization
%%A = S' * diag( data_weights ) * S + diag( regularization_weights );
%%B = S' * diag( data_weights ) * delG_NN + regularization_weights;
%
%% without regularization
%A = S' * diag( data_weights ) * S;
%B = S' * diag( data_weights ) * delG_NN;
%
%new_weights = inv(A) * B;
%
%output_weights( new_weights, data.score_labels( fit_idx ), 0.0 );
%
%% for fitting ref wts and other_scores
%delG_fit =  coeff_matrix * D * new_weights;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%% Fit weights on a subset of score terms only %%%%%%
fit_idx = get_fit_idx( score_terms_to_fit, data.score_labels );
n_weights = length( fit_idx );

%regularization_stdev = 0.02;
regularization_stdev = 0.2;
regularization_weights = ones( n_weights, 1 ) / ( regularization_stdev^2 );
regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'intermol' ) ) ) = 1/0.01^2;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'fa_rep' ) ) )  =  1/0.01^2;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'rna_sugar_close' ) ) )  =  1/0.01^2;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded' ) ) ) = 0.0;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_a' ) ) ) = 0.5;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_c' ) ) ) = 0.5;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_g' ) ) ) = 0.5;
%regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_u' ) ) ) = 0.5;

regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_a' ) ) ) = 0.0;
regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_c' ) ) ) = 0.0;
regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_g' ) ) ) = 0.0;
regularization_weights( find( strcmp( data.score_labels( fit_idx ), 'unfolded_u' ) ) ) = 0.0;

% fitting just fit_idx wts (not other_scores wt)
D = data.scores(:,fit_idx);

data_weights = (1 ./ delG_NN_err.^2);

% try keeping the wt on the "other_scores" fixed
% then the new observable would be delG_NN - other_scores
% other_scores = data.scores(:,1) - sum(data.scores(:,fit_idx),2)
fitting_exp = delG_NN - (data.scores(:,1) - sum(data.scores(:,fit_idx),2));

S = coeff_matrix * D;

%without regularization
A = S' * diag( data_weights ) * S;
B = S' * diag( data_weights ) * fitting_exp;

% with regularization
%A = S' * diag( data_weights ) * S + diag( regularization_weights );
%B = S' * diag( data_weights ) * fitting_exp + regularization_weights;

new_weights = inv(A) * B;

output_weights( new_weights, data.score_labels( fit_idx ), 0.0 );

delG_fit =  coeff_matrix * D * new_weights + (data.scores(:,1) - sum(data.scores(:,fit_idx),2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% Fit only an additional length dependent term %%%%%%
%
%% figure out length of each sequence
%% sequence tags are in seq_labels
%for m = 1:length( seq_labels );
%	len_seq( m ) = length( strfind( char(seq_labels{m}), 'a')) + ...
%	 length( strfind( char(seq_labels{m}), 'c')) + ...
%	 length( strfind( char(seq_labels{m}), 'g')) + ...
%	 length( strfind( char(seq_labels{m}), 'u'));
%end
%
%% try to fit just a length dependent term
%D = len_seq';
%
%data_weights = (1 ./ delG_NN_err.^2);
%
%% for fitting just length dependent term
%fitting_exp = delG_NN - data.scores(:,1);
%%%%%
%
%S = coeff_matrix * D;
%A = S' * diag( data_weights ) * S;
%B = S' * diag( data_weights ) * fitting_exp; % use this for fitting just ref wts, not other_scores wt
%new_weights = inv(A) * B;
%
%new_weights
%
%%output_weights( new_weights, data.score_labels( fit_idx ), 0.0 );
%
%% for fitting just length dependent term:
%delG_fit = coeff_matrix * D * new_weights + data.scores(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots
%%%%%%%%%%%%%%%%%%%%%%%%%%

delG_raw = coeff_matrix * data.scores(:,1);

minE = min( [delG_NN;delG_fit;delG_raw] ) - 1;
%minE = min( [delG_NN;delG_fit] ) - 1;
maxE = max( [delG_NN;delG_fit;delG_raw] ) + 1;
%maxE = max( [delG_NN;delG_fit] ) + 1;

plot( [minE maxE], [minE maxE], 'k-', 'linew',2 );
hold on

refit_color = [1 0.7 0.7];
plot( delG_fit, delG_NN, '.','color',refit_color,'MarkerSize',15);
%plot( delG_fit, delG_NN, '.','color',refit_color );
fprintf( 'fit   rmsd: %8.3f\n', weighted_rmsd( delG_fit, delG_NN, delG_NN_err ) ); 
corrcoef( delG_fit, delG_NN )
corrcoef( delG_raw, delG_NN)

raw_color = [0.7 0.7 1];

%delG_raw = coeff_matrix * data.scores(:,fit_idx);

%plot( delG_raw, delG_NN, '.', 'color', raw_color ); hold on
delG_refit = 0.8*(data.scores(:,1) - sum(data.scores(:,28:31), 2)) + (sum(data.scores(:,28:31), 2)*2.3);
plot( delG_raw, delG_NN, '.', 'color', raw_color,'MarkerSize',15 ); hold on
%plot( 0.5*(delG_raw+15), delG_NN, '.', 'color', raw_color ); hold on
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
