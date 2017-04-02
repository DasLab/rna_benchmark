#!/usr/bin/python

##########################################################
#import imp
#imp.load_source('info_handlers', '../utility/info_handlers.py')
#imp.load_source('file_handlers', '../utility/file_handlers.py')

##########################################################

from sys import exit
from os import system, popen
import subprocess
from os.path import exists, dirname, basename, abspath, isdir
import string
import glob
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
from datetime import datetime
from utility import file_handlers, info_handlers
try:
        from get_sequence import get_sequence, get_sequences_for_res
except Exception, e:
        print e
        print "WARNING ... could not import get_sequences_for_res, required for sequence recovery plotting"
from parse_tag import parse_tag
import re
from difflib import SequenceMatcher


##########################################################

class ScoreData(object):

	__slots__ = ["file", "scores", "score_labels"]

	def __init__(self, file = None):
                self.file = file
		self.scores = None
		self.score_labels = None

##########################################################

class TimeData(object):

	__slots__ = ["inpath", "target", "times", "mean", "stdev"]

	def __init__(self):
		self.inpath = None
		self.target = None
		self.times = []
		self.mean = None
		self.stdev = None

	def times_found(self):
		return ( len(self.times) > 0.0 )

	def update_stats(self):
		if self.times_found():
			self.mean = np.mean( np.array( self.times ) )
			self.stdev = np.std( np.array( self.times ) )
		else:
			self.mean = None
			self.stdev = None

	def get_stats(self):
		return ( self.mean, self.stdev )

	def get_label(self):
		return '%5.0f +/- %4.0f' % self.get_stats()

##########################################################

def get_outfiles( inpath, targets, outfilenames ):
        if not isinstance(targets, list):
                targets = [targets]
	outfiles = []
	for target in targets:
		for outfilename in outfilenames:
			outfiles += glob.glob('/'.join([inpath, target, outfilename]))
	for outfile in outfiles:
		if outfile.find('.out') > 0 and outfile.replace('.out', '.sc') in outfiles:
			outfiles.pop(outfiles.index(outfile))
			continue
		assert( exists( outfile ) )

	return outfiles

##########################################################

def load_score_data( file ):
	data = ScoreData(file=file)
	if not exists( file ):
		return False
	with open( file, 'r' ) as f:
		for l in f:
			if 'SCORE:' not in l: continue
			data.score_labels = l.strip().split()[1:]
			break
	if not data.score_labels:
		return False
	scorenums = [str(x) for x in xrange(2,len(data.score_labels)+2)]
	command = (
		'grep SCORE ' + file
		+ ' | awk \'{print $'+',$'.join(scorenums)+'}\''
		+ ' | grep -v inp'
		+ ' | grep -v R'
		+ ' | grep -v H'
		+ ' | grep -v score'
		+ ' | grep -v pdb'
		+ ' | grep -v descr'
		+ ' | grep -v total'
		#+ ' | grep -v S_'
	)
	data.scores = [s.strip().split() for s in popen(command).xreadlines()]
	return data

##########################################################

def load_data( inpaths, targets, outfilenames ):
	data = {}
	for inpath_idx, inpath in enumerate(inpaths):
		outfiles = get_outfiles( inpath, targets, outfilenames )
		if not len(outfiles):
			continue
                for outfile in outfiles:
                        print 'Reading in ... '+outfile
                score_data_list = map(load_score_data, outfiles)
                score_data_list = filter(None, score_data_list)
                data[ inpath ] = dict([(basename(dirname(score_data.file)), score_data) for score_data in score_data_list])
	return data

##########################################################

def get_date():
	#current_date = popen( 'date' ).read().replace('\n','')
	#current_date = datetime.now().strftime("%b %d, %Y")
	current_date = datetime.now().strftime("%Y-%m-%d")
	return current_date

##########################################################

def get_target_names( target_files, inpaths = None ):
        if inpaths:
                return get_target_names_from_inpaths( inpaths )
	target_names = []
	for file_name in target_files:
		if '..' not in file_name:
                        input_dir = get_path_to_dir(['stepwise_benchmark','benchmark']) + '/input_files/'
                        file_name = input_dir + basename(file_name)
		target_names += get_target_names_from_file( file_name, target_names )
	return target_names

###########################################################

def get_target_names_from_inpaths( inpaths ):
        target_names = []
        for inpath in inpaths:
                target_names += glob.glob('/'.join([inpath, '*']))
        target_names = list(set(map(basename, filter(isdir, target_names))))
        return target_names

###########################################################

def get_target_names_from_file( filename, target_names ):
	assert( exists( filename ) )
        info_fid = file_handlers.TargetDefinitionsFile()
        info_fid.load(open(filename))
        assert( info_fid.validate() )
        return [td.name for td in info_fid.target_definitions]

###########################################################

def get_path_to_dir( dirnames ):
	for dirname in dirnames:
		pwd = popen( 'pwd' ).readline()[:-1].split( '/' )
		for subdir_idx, subdir in enumerate(pwd):
                        if dirname not in subdir:
                                continue
                        return string.join( pwd[:subdir_idx+1], '/' )


###########################################################

def show_times( inpaths, data, target_names, times=None ):
	if not times:
		times_list = get_times( inpaths, data, target_names, verbose=False )
	else:
		times_list = times
	for k in xrange( len(target_names) ):
		if not any ( t[k].times_found() for t in times_list ):
			continue
		print '\n %-31s%8s' % ( target_names[k], 'CPU Time' )
		for n in xrange( len(inpaths) ):
			run_time = times_list[n][k]
			if run_time.times_found():
				print ' Run %d                    %5.0f +/- %4.0f' % ( n+1, run_time.mean, run_time.stdev )
			else:
				print ' Run %d %33s' % ( n+1, 'N/A' )
	print '\n'
	for inpath_idx, inpath in enumerate(inpaths):
		print ' Run %d: %s' % (inpath_idx, basename(inpath))
	return

###########################################################

def get_times( inpaths, data, target_names, verbose=False ):
	time_label = 'time'
	times_list = []
	for inpath_idx, inpath in enumerate(inpaths):
		times = []
		for target in target_names:
			time_data = TimeData()
			time_data.inpath = inpath
			time_data.target = target
			try:
				time_idx = data[ inpath ][ target ].score_labels.index( time_label )
				time_data.times = np.array( [ score[time_idx] for score in data[ inpath ][ target ].scores ], dtype = 'float_' )
			except:
				time_data.times = []
			time_data.update_stats()
			times.append( time_data )
		times_list.append( times )
	if verbose:
		show_times( inpaths, data, target_names, times=times_list )
	return times_list

###########################################################

def get_figure_dimensions( nplots ):
	assert( nplots )
	if nplots <= 5:
		nrows = nplots
	elif nplots <= 9:
		nrows = 3
	elif nplots <= 12:
		nrows = 4
	elif nplots <= 20:
		nrows = 5
	else:
		nrows = 6
	ncols = np.ceil( nplots / float( nrows ) )
	return ( nplots, nrows, ncols )

###########################################################

def setup_figure( nplots, landscape = False ):
	( nplots, nrows, ncols ) = get_figure_dimensions( nplots )
	fig = plt.figure()
	if ( landscape is True ): # landscape
		fig.set_size_inches(11, 8.5)
	else:
		fig.set_size_inches(8.5, 11)
	return ( fig, nplots, nrows, ncols )

###########################################################

def finalize_figure( fig, nplots, nrows, ncols, options=None ):
	# adjust spacing of plots on figure
	if ( nplots == 1 or nrows < ncols ): # landscape
		plt.subplots_adjust(top=.90, bottom=.1,
				    left=.05, right=.98,
				    hspace=.5)
	else:
		plt.subplots_adjust(top=.95, bottom=.1,
				    left=.08, right=.95,
				    wspace=.3, hspace=.5)

	# get date printed to figure
	if options and options.seaborn is True:
                # add legend
                plt.legend(numpoints=1, loc=4, prop={'size':16})
        else:
                # display date in bottom right corner
                plt.figtext(0.95,0.02,
		    get_date(),
		    horizontalalignment='right',
		    fontsize='small')

                # setup global legend based on inpaths
                legend_size = 6
                plot_idx = nplots - 1 if ncols > 1 else nplots
                ax = fig.add_subplot( nrows, ncols, plot_idx)
                legend = ax.legend(bbox_to_anchor=(0., .0, 1., -.225),
                                   loc=9, numpoints=1, prop={'size':legend_size})
       	return

###########################################################

def setup_pdf_page( inpaths, targets, pdfname = None ):
	inpaths = [basename(x) for x in inpaths]
	if not pdfname:
                pdfname = ''
	        #if len( targets ) > 0 and targets[0] != '*': pdfname += string.join(targets, '_') + '_'
	        pdfname += string.join(inpaths, '_vs_') + '.pdf'
                # following was from calebgeniesse -- includes all target names -- huge filename!
		#pdfname = '_'.join(targets)
		#pdfname += '_' + '_vs_'.join(inpaths)
		#pdfname = pdfname if pdfname[0] != '_' else pdfname[1:]
	pdfname += '.pdf' if '.pdf' not in pdfname else ''
	if '/' in pdfname and exists(dirname(pdfname)):
		fullpdfname = pdfname
	else:
		figure_dir = get_path_to_dir(['stepwise_benchmark','benchmark']) + '/Figures/'
		fullpdfname = figure_dir + pdfname
	try:
		pp = PdfPages( fullpdfname )
	except IOError as err:
		if 'File name too long' in err.args:
			datetime_str = datetime.now().strftime("%Y%m%d-%H%M%S")
			fullpdfname = figure_dir + 'make_plots_' + datetime_str + '.pdf'
			pp = PdfPages( fullpdfname )
		else:
			print "\nIOError", err
			exit(-1)
	print '\nMaking figure in: %s\n' % fullpdfname
	return ( pp, fullpdfname )

###########################################################

def get_colorcode( size, options=None ):
        if options and options.seaborn is True:
                return ["#ff0000","#0000ff", "#696969"]
	if size <= 2:
		return [(0.0, 0.0, 0.0, 1.0), (1.0, 0.0, 0.0, 1.0)]
	cmap = plt.get_cmap( 'hot' )
	cnorm = colors.Normalize( vmin=0, vmax=size )
	scalar_map = cmx.ScalarMappable( norm=cnorm, cmap=cmap )
	return [ scalar_map.to_rgba( x ) for x in xrange( size ) ]

###########################################################

def get_title( target ):
	names = get_path_to_dir(['stepwise_benchmark','benchmark']) + '/scripts/python/plot_util/titles.txt'
	try:
		lines = open( names, 'r' ).readlines()
	except:
		return target
	for line in lines:
		if target in line:
			return line.replace('\n','').split(': ')[1]
	return target

###########################################################

def get_fasta_from_silent_file(outfile):
        """ Returns .fasta file listing sequencename and base sequences provided from outfile """

        fasta_file = outfile+".fasta"
        fasta = open(fasta_file,'w')

        tags, sequences, n_res = get_sequences_from_silent_file(outfile)

        for tag, sequence in zip(tags,sequences):
                fasta.write(">{}\n{}\n".format(tag,sequence))
        fasta.close()
        return fasta_file, n_res

def get_sequences_from_silent_file(outfile,tags=None):

        original_sequence = ""
        seqname = ""
        res_num, res_chain = None, None
        index = []
        sequences = []
        good_tags = []
        other_pose = False
        n_res = None

        for line in open(outfile).readlines():
                line = line.strip()

                if "FULL_SEQUENCE" in line:

                        line = line.split()
                        full_seq_index = line.index("FULL_SEQUENCE")
                        original_sequence = line[full_seq_index+1]

                        if "CONVENTIONAL_RES_CHAIN" in line:
                                res_chain_index = line.index("CONVENTIONAL_RES_CHAIN")
                                conventional_res_chain = line[res_chain_index+1]
                                conventional_res, conventional_chain = parse_tag(conventional_res_chain)
                                #print "conventional res chain:", zip(conventional_res,conventional_chain)

                        index_list = [n_char.start() for n_char in re.finditer('n', original_sequence)]

                        # show neighboring base pairs of design sequence
                        index_list += [idx-1 for idx in index_list if idx-1 not in index_list]
                        index_list += [idx+1 for idx in index_list if idx+1 not in index_list]
                        index_list.sort()

                        n_res = [conventional_res[i] for i in index_list]
                        continue

                if "ANNOTATED_SEQUENCE:" in line:

                        line = line.split()
                        if tags and seqname not in tags:
                                continue

                        if seqname != line[-1] or other_pose:
                                continue

                        sequence = re.sub(r'\[.*?\]', '', line[1])

                        if len(original_sequence) != len(sequence):

                                filled_sequence = ""
                                for res, char in zip(conventional_res, original_sequence):
                                        if res not in res_num:
                                                filled_sequence += 'n'
                                        else:
                                                filled_sequence += sequence[res_num.index(res)]

                                sequence = filled_sequence

                        assert (len(sequence) == len(original_sequence))

                        sequence = "".join([sequence[i] for i in range(len(sequence)) if i in index_list])

                        sequences.append(sequence)
                        good_tags.append(seqname)
                        continue

                if "RES_NUM" in line:

                        line = line.split()

                        other_pose = (seqname == line[-1])
                        if not other_pose:
                                seqname = line[-1]
                                if tags and seqname not in tags:
                                        continue
                                res_input = ",".join(line[1:-1])
                                res_num, res_chain = parse_tag(res_input)

        return good_tags,sequences, n_res

###################################################


def get_sequence_recovery(inpath, target, outfilenames, tag):

        outfilenames = [name for name in outfilenames if not name.endswith(".sc")]
        outfiles = get_outfiles(inpath, target, outfilenames)

        if not len(outfiles):
                return 0

        tags, sequences, n_res = get_sequences_from_silent_file(outfiles[0],tags = [tag])

        if len(sequences) == 0:
                return 0

        sequence = sequences[0]
        target_pdbs = glob.glob('/'.join([inpath, target, target+"*RNA.pdb"]))
        native_pdb = [ pdb for pdb in target_pdbs if "START" not in pdb and "HELIX" not in pdb]
        if not len(native_pdb):
                return 0
        n_res = " ".join(map(str,n_res))
        native = "".join(get_sequences_for_res(native_pdb[0],n_res))

        seq_recovery = 1
        # go through letters of native sequence
        # mark 1 if match
        # len(sequence)
        # match/sequence * 100
        seq_recovery = SequenceMatcher(None, sequence, native).ratio()*100
        print sequence
        print native, seq_recovery
        return seq_recovery
