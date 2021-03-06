#!/usr/bin/python

###############################################################################
### imports
###############################################################################
import sys
import os
from os.path import exists, dirname, basename, abspath, isdir
import stat
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from make_plots_util import *
import subprocess
from matplotlib.font_manager import FontProperties
import argparse
import re
from parse_tag import parse_tag

###############################################################################
### file functions
###############################################################################

def generate_weblogo(inpath, target, outfile):
        """ Returns .weblogo.png image of sequence given by outfile """

        outfile = "/".join([inpath,target,outfile])
        fasta_file, n_res = get_fasta_from_silent_file(outfile)  
        res_csv = ','.join(map(str,n_res))
       
        if os.stat(fasta_file).st_size == 0:
                return False

        weblogo_file = outfile.replace(".out",".weblogo.png")

        cmd = [
        "weblogo", 
        "-F png",
        "--resolution 600",
        "--alphabet 'CGAUN'",
        "--annotate '{}'".format(res_csv),
        "--size large",
        "--color green cC 'Cytosine'",
        "--color red gG 'Guanine'",
        "--color orange aA 'Adenine'",
        "--color blue uU 'Uracil'",
        "--color black nN 'N'",
        "--errorbars NO",
        "< {} > {}".format(fasta_file, weblogo_file)]

        cmd = " ".join(cmd)
        os.system(cmd)

        return weblogo_file

###############################################################################
### main function
###############################################################################
def make_weblogos(argv):

	# initialize options
	options = init_options(argv)
	inpaths = options.inpaths
	silentfile = options.silentfile
	target_files = options.target_files
	targets = options.targets
        pdfname = options.pdfname

	# check options
	inpaths = [abspath(x) for x in inpaths if exists(x) and isdir(x)]
        assert(len(inpaths) == 1)

	inpath = inpaths[0]

        if targets[0] != '*':
                targets = targets
        elif target_files is not None:
                targets = get_target_names( target_files )
        else:
                targets = get_target_names( target_files, inpaths )
	for target in targets:
		print "Target: "+target
                
        nplots = len(targets)
        if nplots >= 3:
                nplots = 3
        
        if pdfname is None:
                pdfname = basename(inpath)+'_weblogos.pdf'

        # setup pdf and figure, return handles
	( pp, fullpdfname ) = setup_pdf_page( inpaths, targets, pdfname )
	( fig, nplots, nrows, ncols ) = setup_figure( nplots)
        
	# iterate over targets
        count = 1
	for plot_idx, target in enumerate(targets, start=1):
                
                if count > nplots or count == 1 or plot_idx == len(targets):
                        if count != 1:
                                count = 1
                                # get date printed to figure
                                plt.figtext(0.95,0.98,
                                            basename(inpath),
                                            horizontalalignment='right',
                                            fontsize='large')
                                pp.savefig(DPI=1200)
                        ( fig, nplots, nrows, ncols ) = setup_figure( nplots )
                        
                weblogo = generate_weblogo(inpath, target, silentfile)
                #print weblogo
                
                if weblogo is False:
                        continue
                
                # plot image on subplot
                ax = fig.add_subplot( nrows, ncols, count )
                image = mpimg.imread(weblogo)

                plt.imshow(image)
                ax.set_title( target, fontsize=8, weight='bold' )
                
                count += 1
                plt.axis('off')

	# finalize (adjust spacing, print date)
	#finalize_figure( fig, nplots, nrows, ncols )

	# save as pdf and close
	pp.savefig(DPI=1200)
	pp.close()

	# open pdf
	out, err = subprocess.Popen(['uname'], stdout=subprocess.PIPE).communicate()
	if 'Darwin' in out:
		subprocess.call(['open',fullpdfname])
	if 'Linux' in out:
		subprocess.call(['xdg-open',fullpdfname])

	return True

###############################################################################
### initialization functions
###############################################################################
def init_options(argv):
    if isinstance(argv, argparse.ArgumentParser):
    	return argv
    options = init_options_parser().parse_args(args = argv)
    return options

def init_options_parser():
	parser = argparse.ArgumentParser(description='Plot scores from silent files.')
	parser.add_argument(
		'inpaths',
		nargs='+',
		help='List of paths too silent files.'
	)
	parser.add_argument(
		'-silentfile',
	    help='Name of silent file.',
	    default='swm_rebuild.out'
	)
	parser.add_argument(
		'-target_files',
	    nargs='+',
	    help='List of additional target files.',
	    default=None
	)
	parser.add_argument(
		'-targets',
	    nargs='+',
	    help='List of targets.',
	    default=['*']
	)
	
	parser.add_argument(
		'-o','--pdfname',
	    help='File name to save as pdf.',
	    default=None
	)
	return parser


###############################################################################
### main
###############################################################################
if __name__=='__main__':

	sys.exit(make_weblogos(sys.argv[1:]))

