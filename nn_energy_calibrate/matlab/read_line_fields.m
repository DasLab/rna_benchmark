function line_fields = read_line_fields( filename, delimiter );
%  line_fields = read_line_fields( filename );

if ~exist( 'delimiter','var') delimiter = '\t'; end;
fid = fopen( filename );
count = 0;
while ~feof( fid );
  line = fgetl( fid );
  lines = strsplit( line, '\r' );
  for m = 1:length( lines )
    line = lines{m};
    if length( strsplit( line ) ) == 1; continue; end;
    %if length( strsplit( line ) ) == 0; continue; end;
    if line(1) == '%'; continue; end
    if line(1) == '#'; continue; end
    line  = strrep( line, '"','');
    count = count+1;
    line_fields{count} = strsplit( line, delimiter );
  end
end
