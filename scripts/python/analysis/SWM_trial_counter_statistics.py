#!/usr/bin/python

###############################################################################

from glob import glob
import subprocess as sp
from os.path import isdir, basename, exists, abspath
import os
import subprocess
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

###############################################################################

FIGURES_DIR = abspath(__file__).split('benchmark/')[0]+'benchmark/Figures/'

###############################################################################

def extract_trial_info( line ):
    cols = [c.strip() for c in line.split(':')[-1].split(';')]
    try:
        info = [float(v.split('=')[-1].split()[0]) for v in cols]
        return ( int(info[0]), info[1], info[2] )
    except:
        print cols
    return (None, None, None)

###############################################################################

def get_trial_counter_statistics( outfile ):

    dirs = glob('swm_legacy_match_add_pdf*')+['swm_legacy_match']
    dirs.sort(key=len)
    targets = map(basename, filter(isdir, glob('swm_legacy_match/*')))
    targets.sort()
    trials_info = {}

    for dir in dirs:
        print "Checking logfiles in: ", dir
        trials_info[ dir ] = {}
        for target in targets:
            print "Checking logfiles in: ", dir+'/'+target
            trials_info[dir][target] = [[],[],[]]
            lines = []
            logfiles = glob('%s/%s/[0-9]*.out' % (dir, target))
            for logfile in logfiles:
                with open( logfile, 'r' ) as f:
                    for line in f:
                        if not 'TrialCounter:' in line: continue
                        if not 'add trials=' in line: continue
                        lines.append(line)
            for idx, line in enumerate(lines):
                if idx == 0: continue
                if 'NO ACCEPTS.' in line: continue
                if extract_trial_info( line )[0] != 1: continue
                trial_info = extract_trial_info( lines[idx-1] )
                if not trial_info[0] or trial_info[0] == 1: continue
                for idx, info in enumerate(trial_info):
                    trials_info[dir][target][idx].append( info )

    fid = open( outfile, 'w' )
    fid.write('TARGET\tINPATH\tPDF\tTRIALS\tACCEPTS\tEDROP_PER_TRIAL\n')
    for target in targets:
        print "\nTARGET:",target
        for dir in dirs:
            pdf = int(dir.split('_')[-1]) if 'pdf' in dir else 1
            aves = [sum(info)/len(info) for info in trials_info[dir][target]]
            print "INPATH: %-40s TRIALS: %-4d ACCEPTS: %-8f ENERGY_DROP/TRIAL: %-8f"%(
                dir, aves[0], aves[1], aves[2])
            fid.write("%s\t%s\t%d\t%d\t%f\t%f\n"%(
                    target,dir,pdf,aves[0],aves[1],aves[2]))
    fid.close()

    return

###############################################################################

def plot_trial_counter_statistics( outfile ):

    # get data
    pdf_data = []
    acc_data = []
    edrop_data = []

    if not exists( outfile ):
        get_trial_counter_statistics( outfile )
    assert( exists( outfile ) )

    with open( outfile, 'r' ) as f:
        for line in f:
            cols = line.split()
            if 'TARGET' in line:
                pdf_idx = cols.index('PDF')
                acc_idx = cols.index('ACCEPTS')
                edrop_idx = cols.index('EDROP_PER_TRIAL')
                continue
            pdf_data.append( int(cols[pdf_idx]) )
            acc_data.append( float(cols[acc_idx]) )
            edrop_data.append( float(cols[edrop_idx]) )

    # setup pdf
    fullpdfname = outfile.replace('.out','.pdf')
    pp = PdfPages( fullpdfname )
    print '\nMaking figure in', fullpdfname

    # setup figue
    fig = plt.figure(1)
    titlesize = 14
    labelsize = 12
    ticklabelsize = 10

    # plot pdf vs acceptance u 3:5
    ax = fig.add_subplot(2, 1, 1)
    ax.plot( pdf_data,
             acc_data,
             marker='o',
             markersize=4,
             linestyle=' ')
    ax.set_title('SWM Add Proposal Density Factor vs. Add Move Acceptance',
                 fontsize=titlesize
                 )
    ax.set_xlabel('proposal density factor', fontsize=labelsize)
    ax.set_ylabel('acceptance', fontsize=labelsize)
    ax.set_xscale('log')
    ax.set_xlim(1,100000)
    for tl in ax.yaxis.get_ticklabels()+ax.xaxis.get_ticklabels():
        tl.set_fontsize(ticklabelsize)

    # plot pdf vs edrop/trial u 3:6
    ax = fig.add_subplot(2, 1, 2)
    ax.plot( pdf_data,
             edrop_data,
             marker='o',
             markersize=4,
             linestyle=' ')
    ax.set_title('SWM Add Proposal Density Factor VS. Add Move Energy Drop/Trial',
                 fontsize=titlesize
                 )
    ax.set_xlabel('proposal density factor', fontsize=labelsize)
    ax.set_ylabel('energy drop/trial', fontsize=labelsize)
    ax.set_xscale('log')
    ax.set_xlim(1,100000)
    for tl in ax.yaxis.get_ticklabels()+ax.xaxis.get_ticklabels():
        tl.set_fontsize(ticklabelsize)

    # finalize figure
    plt.subplots_adjust(bottom=.1,left=.1, right=.95, top=.9, wspace=.3, hspace=.5)

    # close fig / open pdf
    pp.savefig()
    pp.close()

    uname = os.uname()
    if 'Darwin' in uname:
        subprocess.call(['open',fullpdfname])
    if 'Linux' in uname:
        subprocess.call(['xdg-open',fullpdfname])

    return


###############################################################################

if __name__=='__main__':

    import argparse

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p','--plot',action='store_true')
    parser.add_argument('-o','--outfile',default="TrialCounterStats.out")
    args=parser.parse_args()

    if args.plot:
        plot_trial_counter_statistics( args.outfile )
    else:
        get_trial_counter_statistics( args.outfile )
