function [ tags_NN, delG_NN, delG_NN_err, coeff_matrix ] = get_NN_linear_combinations( linear_combination_file, tags );

for i = 1:length( tags )
  %tags{i} = basename( tags{i} );
  tags{i} = strrep( tags{i}, '_helix','' );
  tags{i} = strrep( tags{i}, '.pdb','' );
end

all_line_fields = read_line_fields( linear_combination_file );


count = 0;
for n = 1:length( all_line_fields )

  line_fields = all_line_fields{n};

  if length( line_fields ) < 1; continue; end
  if length( line_fields{1} ) < 1; continue; end
  
  count = count + 1;
  
  tags_NN{count}     = line_fields{1};
  delG_NN(count)     = str2num( line_fields{2} );
  delG_NN_err(count) = str2num( line_fields{3} );
  coeff_matrix( count, :) =  zeros( 1, length(tags) );

  for m = 4 : 2 : length( line_fields )-1
  %for m = 4 : 2 : length( line_fields )    
    coeff = str2num( line_fields{m} );
    
    tag = line_fields{m+1} ;

    % 04-11-16: I don't know why this isn't working... just do it the long way instead (matlab v2014a)
    %idx = find(strcmp( tags, tag));
    %if ( length(idx) ~= 1 ) fprintf( '%s --> %d\n',  line_fields{m+1}, length(idx) ); error( 'tag problem' ); end;
    idx = 0;
    i = 0;
    for k = 1:length(tags)
      i = i + 1;
      if strcmp(tags{k},tag); idx = i; end;
    end;

    if ( idx == 0 ) fprintf( '%s --> %d\n',  line_fields{m+1}, length(idx) ); error( 'tag problem' ); end;

    coeff_matrix( count, idx ) = coeff_matrix( count, idx ) + coeff;  
  end

end

delG_NN = delG_NN';
