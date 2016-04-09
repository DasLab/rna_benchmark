function data = add_separate_unfolded_energies( data, tags, scorefilename ); 

rna_char = 'gacu';

ref_weights = get_ref_weights( scorefilename );

for n = 1:length(rna_char)
  r = rna_char(n);
  unfolded_score_labels{n} = ['unfolded_',r];
  for m = 1:length( tags );
    unfolded( m, n ) = ref_weights(n) * length( strfind( tags{m}, r ) );
  end
end

data.scores = [ data.scores, unfolded ];
data.score_labels = [ data.score_labels, unfolded_score_labels ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ref_weights = get_ref_weights( scorefilename );

all_line_fields = read_line_fields( scorefilename, ' ' );
for n = 1:length( all_line_fields)
  line_fields = all_line_fields{n};  
  if length( line_fields ) > 1 && strcmp( line_fields{1}, 'METHOD_WEIGHTS' );
    method_weights = str2num( char(line_fields( end-3: end )) );
  end
  if length( line_fields ) > 1 && strcmp( line_fields{1}, 'ref' );
    ref_weight = str2num( line_fields{2} );
  end
end

ref_weights = ref_weight * method_weights;
%ref_weights = 0.95 * [2.64 2.28 1.8 2.4];
