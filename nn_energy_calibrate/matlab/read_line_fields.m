function line_fields = read_line_fields( filename, delimiter );
%  line_fields = read_line_fields( filename );

if ~exist( 'delimiter','var') delimiter = '\t'; end;
fid = fopen( filename );
count = 0;
while ~feof( fid );
  line = fgetl( fid );
  lines = split_string( line, '\r' );
  for m = 1:length( lines )
    line = lines{m};
    if length( split_string( line ) ) == 0; continue; end;
    if line(1) == '%'; continue; end
    if line(1) == '#'; continue; end
    line  = strrep( line, '"','');
    count = count+1;
    line_fields{count} = split_string( line, delimiter );
  end
end
