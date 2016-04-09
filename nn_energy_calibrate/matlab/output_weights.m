function output_weights( X, feature_names, cutoff, tag )

if isempty( X ); return; end;
if ~exist( 'tag', 'var' ) tag = 'WEIGHTS'; end;

fprintf( '\n')
fprintf( [tag,'\n'])

for i = 1:length(X);    
  if abs( X(i) ) >= cutoff
    fprintf( '%8.2f  %s\n', X(i), feature_names{i} ) ;   
  end
end
fprintf( '\n')
