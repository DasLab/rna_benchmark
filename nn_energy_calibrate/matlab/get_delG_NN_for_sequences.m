function delG = get_delG_NN_for_sequences( seq1, seq2, turner_rules );

rna_char = turner_rules.rna_char;

delG = turner_rules.init;

assert( length( seq1 ) == length( seq2 ) );
seq1_helix = seq1;
seq2_helix = seq2;

% take care of termini
if is_blank_seq( seq1(1) ) | is_blank_seq( seq2(end) )
  seq1_helix = seq1_helix(2 : end);
  seq2_helix = seq2_helix(1 : (end-1) );
  if ( is_blank_seq( seq1(1) ) )
    num1 = strfind( rna_char, seq2( end-1 ) );
    num2 = strfind( rna_char, seq2( end ) );
    delG = delG + turner_rules.dangle_3prime( num1, num2 );
  else 
    assert( is_blank_seq( seq2(end) ) );
    num1 = strfind( rna_char, seq1( 2 ) );
    num2 = strfind( rna_char, seq1( 1 ) ); % that's the table orientation in Serra & Turner
    delG = delG + turner_rules.dangle_5prime( num1, num2 );    
  end
end

if is_blank_seq( seq1(end) ) | is_blank_seq( seq2(1) )
  seq1_helix = seq1_helix(1 : (end-1) );
  seq2_helix = seq2_helix(2 : end );
  if ( is_blank_seq( seq1(end) ) )
    num1 = strfind( rna_char, seq2( 2 ) );
    num2 = strfind( rna_char, seq2( 1 ) );
    delG = delG + turner_rules.dangle_5prime( num1, num2 );
  else 
    assert( is_blank_seq( seq2(1) ) );
    num1 = strfind( rna_char, seq1( end-1 ) );
    num2 = strfind( rna_char, seq1( end ) );
    delG = delG + turner_rules.dangle_3prime( num1, num2 );    
  end
end

% nearest neighbor rules.
rna_bps = turner_rules.rna_bps;
i = 1;
while ( i < length( seq1_helix ) )

  if ( i+3 <= length( seq1_helix ) & ...
       strcmp( seq1_helix(i:i+3), 'gguc' ) & ...
       strcmp( seq2_helix( end-i-2 : end-i+1), 'gguc' ) )
    delG = delG + turner_rules.special_gguc;
    i = i+3;
    continue
  end

  num1 = find( strcmp( rna_bps, [seq1_helix(i),   seq2_helix( end-i+1) ]  ) );
  num2 = find( strcmp( rna_bps, [seq1_helix(i+1), seq2_helix( end-i  ) ] ) );
  delG = delG + turner_rules.delG_NN( num1, num2 );
  i = i+1;
end


for i = [1 length(seq1_helix) ];
  terminal_bp  = [seq1_helix(i), seq2_helix(end-i+1)];
  if ~isempty( find( strcmp(  terminal_bp, {'au','ua'} ) ) )
    delG = delG + turner_rules.per_term_AU;
  end
end

