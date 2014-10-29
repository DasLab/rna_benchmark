#!/usr/bin/python

##########################################################

from os.path import exists, dirname, basename
from os import popen
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from make_plots_util import load_score_data, get_target_names, jet

##########################################################

def make_plots( inpaths, outfilename='swm_rebuild.out', colorcode=None ):
  
  nrows = 4
  ncols = 3
  if not colorcode: colorcode = [ (0.0, 0.0, 0.0, 1.0), (1.0, 0.0, 0.0, 1.0) ]
  if len( colorcode ) < len( inpaths ): colorcode = jet( len( inpaths ) )

  target_names = get_target_names()
  for target_name in target_names:
    print target_name
  
  data = []
  which_target = []
  outfiles_list = []
  for inpath in inpaths:
   
    assert( exists( inpath ) )
    outfiles = popen( 'ls -1 '+inpath+'/*/'+outfilename ).read().split('\n')[:-1] 
   
    for outfile in outfiles:  
      print 'Reading in ... '+outfile
      assert( exists( outfile ) )

    which_target.append( map( lambda x: target_names.index( basename( dirname( x ) ) ), outfiles ) )
    data.append( map( lambda x: load_score_data( x ), outfiles ) )
    outfiles_list.append( outfiles )


  '''
  ###################################################

  time_name = 'time'
  Times = []
  
  for n in xrange( len( inpaths ) ):  
    times = []
    for k in xrange( len( outfiles[n] ) ):
      if not len(data[n][k].scores): continue
      time_idx = data[n][k].score_labels.index( time_name )
      times.append( [ score[time_idx] for score in data[n][k].scores ] ) 
      
    Times.append( times )  
  
  #print Times
  #print len(Times)
  #for times in Times:
  #  print len(times)

  print 
  print '%30s' % 'target'
  for n in xrange( len( inpaths ) ):
    print 'Run %d' % n
    print  
 
    k = 0
    for outfile in outfiles_list[ n ]:
      k+=1
      print '%30s' % target_names[ which_target[n][k] ] 
    
    for n in xrange( len( inpaths ) ):
      mean_time = 0#np.mean( Times[n][k] )
      std_time = 0#np.std( Times[n][k] )
      print '   %5.0f +/- %4.0f' % ( mean_time, std_time )
    print
  print
  for n in xrange( len( inpaths ) ):
    print ' Run %d: %s' % ( n, inpaths[n] )
  
  '''
  
  ###################################################
  
  pdfname = basename( inpaths[0] )
  if len( inpaths ) > 1:
    for k in xrange( 1, len( inpaths ) ): pdfname += '_vs_' + basename( inpaths[k] )
  fullpdfname = 'Figures/' + pdfname + '.pdf'
  print '\nMaking figure in: %s\n' % fullpdfname 
  pp = PdfPages( fullpdfname )
  
  plt.figure(1)
  titles = []
  
  for n in xrange( len( inpaths ) ):
    score_name = 'score'
    rms_name = 'rms_fill'
      
    print which_target[n]
    print len( outfiles_list[n] )

    for k in xrange( len( outfiles_list[n] ) ):

      plt.subplot( nrows, ncols, np.mod( which_target[n][k] -1, nrows*ncols ) + 1 )
      if not len( data[n][k].scores ): continue 
      
      score_idx = data[n][k].score_labels.index( score_name )
      rms_idx = data[n][k].score_labels.index( rms_name )
      score_data = [ score[score_idx] for score in data[n][k].scores ]  
      rms_data = [ score[rms_idx] for score in data[n][k].scores ]  

      plt.plot( rms_data, score_data, marker='.', markersize=3, color=colorcode[n], linestyle=' ' )
      plt.title( target_names[ which_target[n][k] ] )
      plt.xlim( 0, 12 )

      if ( np.mod( which_target[n][k], ncols ) == 1 ):  plt.ylabel( score_name )
      if ( np.floor( (which_target[n][k]-1) / ncols ) == nrows-1 ): plt.xlabel( rms_name )
      #if n == len( inpaths ):
        ################
        #ylim0 = plt.ylim()
        #line( [ 1 1 ], ylim0,'color','k','linestyle',':','selectionhighlight','off')
        #line( [ 2 2 ], ylim0,'color','k','selectionhighlight','off')
        #plt.ylim( ylim0 );
        ############

    titles.append( basename( inpaths[n] ) )
  
  plt.subplot(nrows,ncols,1)
  plt.legend( titles, prop={'size':6} )
 
  pp.savefig()
  pp.close()
  plt.show()

  return

##########################################################

if __name__=='__main__':

  import argparse

  parser = argparse.ArgumentParser(description='Make plots of scores from silent files.')
  parser.add_argument('inpaths', nargs='+', help='List of paths to silent files.')
  parser.add_argument('-outfilename', help='Name of silent file.', default='swm_rebuild.out')
  args=parser.parse_args()

  make_plots( args.inpaths, outfilename=args.outfilename )