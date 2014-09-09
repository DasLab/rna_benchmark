function d = make_plots( inpaths, colorcode);

nrows = 4; ncols = 3;
if ~exist( 'colorcode', 'var' ) colorcode = [0 0 0; 1 0 0]; end;
if size( colorcode, 1 ) < length( inpaths ); colorcode = jet( length( inpaths ) ); end;
clf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target_names = get_target_names();
for n = 1:length( inpaths )
  inpath = inpaths{n};
  assert( exist( inpath,'dir' )>0 );

  outfilename = 'swm_rebuild.out';
  outfiles = split_string( ls( '-1', [inpath,'/*/',outfilename ] ), '\n' );
  for  k= 1:length( outfiles )
    fprintf( ['Reading in... ', outfiles{k}, '\n'] );
    dirn = dirname( outfiles{k} );
    target = basename( dirn(1:end-1) );
    which_target{n,k} = find( strcmp( target_names, target ) );
    [data{n,k}, tags{n,k} ] = load_score_data( outfiles{k} );
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_name = 'time';
for n = 1:length( inpaths )
  for  k= 1:length( outfiles )
    %subplot( nrows, ncols, mod( which_target{n,k} -1, nrows*ncols ) + 1 );
    if length( data{n,k} ) == 0; continue; end; 
    time_idx = find(strcmp( data{n,k}.score_labels, time_name ));
    %T = [0:50:20000];
    %h = hist( data{n,k}.scores(:,time_idx), T );
    %plot( T, h, 'color',colorcode(n,:) );
    %hold on;
    %title( target_names{ which_target{n,k} },'interp','none' );
    %xlim( [0 20000] );
    %if ( mod( which_target{n,k}, ncols ) == 1 ); ylabel( 'N' ); end;
    %if ( floor( (which_target{n,k}-1)/ncols) == nrows-1 ); xlabel( time_name,'interp','none' ); end;
    times{n,k} = data{n,k}.scores(:,time_idx );
  end  
end

fprintf( 1, '\n' )
fprintf( 1, '%30s', 'target' );
for n = 1:length( inpaths )
  fprintf( 1, '            Run %d', n )
end
fprintf( 1, '\n' );

for  k= 1:length( outfiles )
  fprintf( 1, '%30s', target_names{ which_target{n,k} } );
  for n = 1:length( inpaths )
    mean_time = mean( times{n,k} );
    std_time = std( times{n,k} );
    fprintf( 1, '   %5.0f +/- %4.0f', mean_time, std_time )
  end
  fprintf( 1, '\n' );
end
fprintf( 1, '\n' );
for n = 1:length( inpaths )
  fprintf( ' Run %d: %s\n', n, inpaths{n} );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(figure(1), 'position', [300 200 500 600] );
for n = 1:length( inpaths )
  score_name = 'score'; rms_name = 'rms_fill';
  for  k= 1:length( outfiles )
    subplot( nrows, ncols, mod( which_target{n,k} -1, nrows*ncols ) + 1 );
    if length( data{n,k} ) == 0; continue; end; 
    score_idx = find(strcmp( data{n,k}.score_labels, score_name ));
    rms_idx = find(strcmp( data{n,k}.score_labels, rms_name ));
    plot( data{n,k}.scores(:,rms_idx), data{n,k}.scores(:,score_idx),'.','color',colorcode(n,:) );
    hold on;
    title( target_names{ which_target{n,k} },'interp','none' );
    xlim( [0 12] );
    if ( mod( which_target{n,k}, ncols ) == 1 ); ylabel( score_name ); end;
    if ( floor( (which_target{n,k}-1)/ncols) == nrows-1 ); xlabel( rms_name,'interp','none' ); end;
    if n == length( inpaths )
      ylim0 = ylim();
      line( [ 1 1 ], ylim0,'color','k','linestyle',':','selectionhighlight','off')
      line( [ 2 2 ], ylim0,'color','k','selectionhighlight','off')
      ylim( ylim0 );
    end
  end
  titles{n} = basename( inpaths{n} );
end

set(gcf, 'PaperPositionMode','auto','color','white');

subplot(nrows,ncols,1); 
h = legend( titles,1 );
set(h ,'interp','none','fontsize',6)

if length ( inpaths ) > 1
  pdfname = basename( inpaths{1} );
  for k = 2:length( inpaths );  pdfname = [ pdfname, '_vs_',basename( inpaths{k} ) ];  end
  fullpdfname = ['Figures/',pdfname, '.pdf'];
  fprintf( '\nMaking figure in: %s\n\n', fullpdfname );
  export_fig( fullpdfname );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function target_names = get_target_names( );

target_names = get_target_names_from_file( '../favorites.txt', {} );
target_names = get_target_names_from_file( '../favorites2.txt', target_names );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function target_names = get_target_names_from_file( filename, target_names );

fid = fopen( filename );
line = fgetl( fid );
while ~feof( fid )
  line = fgetl( fid );
  cols = split_string( line );
  if length( cols ) == 0; continue;end;
  target_names = [ target_names, cols{1} ];
end
fclose( fid );
