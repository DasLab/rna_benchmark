function [ seq_labels, delG_NN ] = get_delG_NN_for_tags( tags, turner_rules );

seq_labels = {}; delG_NN = [];

delG_NN = zeros( length(tags ), 1 );
for i = 1:length( tags )
  pdb_file = basename( tags{i} );

  cols = split_string( pdb_file, '_' );
  seq1 = cols{1};
  seq2 = cols{2};
  
  delG_NN(i) = get_delG_NN_for_sequences( seq1, seq2, turner_rules );
  seq_labels{i} = [seq1,'/',seq2];

end

