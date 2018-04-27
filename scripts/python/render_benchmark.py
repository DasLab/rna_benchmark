#!/usr/bin/python


###############################################################################
### imports
###############################################################################
import sys
import os
from os.path import exists, dirname, basename, abspath, isdir
import stat
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'

import matplotlib.image as mpimg
import PIL.Image as Image
import matplotlib.gridspec as gridspec
import numpy as np
from make_plots_util import *
import subprocess
from matplotlib.font_manager import FontProperties
import argparse
import re
from parse_tag import parse_tag
from glob import glob

import inspect
__source_file__ = os.path.splitext(os.path.basename(
    inspect.getsourcefile(inspect.currentframe())))[0]
__source_dir__ = os.path.dirname(
    inspect.getsourcefile(inspect.currentframe()))


###############################################################################
### globals
###############################################################################
#plt.rc('text', usetex=True)
plt.rc('text', usetex=False)
print 'usetex:', plt.rcParams['text.usetex']



###############################################################################
### helpers
###############################################################################
def format_tex(*s, **kw):
    delim = ' ' if 'delim' not in kw else kw['delim']
    encode = True if 'encode' not in kw else kw['encode']
    bold = False if 'bold' not in kw else kw['bold']
    italics = False if 'italics' not in kw else kw['italics']
    if type(s) != str:
        s = delim.join(s)
    if plt.rcParams['text.usetex'] is False:
        return s
    if bold:
        s = '\\textbf{' + s + '}'
    if italics:
        s = '\\textit{' + s + '}'
    if encode is False:
        return s
    s = s.replace('_','\_')
    #s = s.replace('-','\-')
    #s = s.replace('/','\/')
    return r"%s" % s

def raw_string(s):
    if isinstance(s, str):
        s = s.encode('string-escape')
    elif isinstance(s, unicode):
        s = s.encode('unicode-escape')
    return s

def fixup_description(desc, props):
	s = ''.join([
		' (',
		format_tex('strands:', props['strands'], bold=False, encode=False),
		', ',
		format_tex('PDB:', props['pdb'], bold=False, encode=False),
		')'
	])
	desc += format_tex(s, bold=False, italics=False, encode=False)
	return desc

def parse_descriptions(desc_file):
	if desc_file is None:
		return None
	with open(desc_file, 'r') as f:
		descriptions = dict([tuple(l.strip().split('\t')[:2]) for l in f.readlines()[1:]])
		for k, v in descriptions.iteritems():
			if v != '-':
				continue
			descriptions[k] = v.replace('-', '')

	props_file = desc_file.replace('descriptions', 'motif-props')
	motif_props = None
	if os.path.exists(props_file):
		with open(props_file, 'r') as f:
			lines = [l.strip().split('\t') for l in f.readlines()[1:] if len(l.strip()) and "Name" not in l]
			motif_props = [
						   tuple([
								  l[0], 
								  dict([('pdb',l[1]), ('strands',l[2])])
								  ]) for l in lines if len(l) > 2
						   ]
			motif_props = dict(motif_props)
		
		# append motif props
		for name, desc in descriptions.iteritems():
			if name not in motif_props:
				continue
			mp = motif_props[name]
			if mp['pdb'] == '--':
				mp['pdb'] = 'n/a'
			descriptions[name] = fixup_description(desc, mp)
				
	for k, v in descriptions.iteritems():
		if 'fixed' not in k:
			continue
			descriptions[k] += ', modeled with rigid-body constraints'
	
	
	# for testing
	for k, v in descriptions.iteritems():
		if len(v):
			continue
		desc = ''.join(['X' for x in xrange(40)])#"This is an example description."
		desc = '\n'.join([desc for x in range(0,2)])
		descriptions[k] = desc
	
	#for k,v in descriptions.iteritems():
		#print k
		#print v
		#print
	
	return descriptions

def add_description(fig, ax_desc, desc, text_obj=None, char_w=None):

	if char_w is None:
		char_w = max([len(s) for s in desc.split('\n')])
	else:
		desc_full = desc.replace('\n',' ')
		desc = ''
		last_space = 0
		split_ok = True
		for idx, x in enumerate(desc_full):
			if x in ['(', '{']:
				split_ok = False
			if x in [')', '}']:
				split_ok = True
			if x == ' ' and split_ok:
				last_space = idx
			if idx > 0 and idx % char_w == 0:
				if x == ' ' and split_ok:
					desc += '\n'
					continue
				elif last_space < len(desc):#else:
					desc = bytearray(desc)
					desc[last_space] = '\n'
					desc = str(desc)
					#if desc.count('\n') < 2:
					#    desc_w = len(desc[:last_space])
					#    nspace = int(np.floor((char_w - desc_w)/2))
					#    desc = ''.join([' ' for _ in range(0, nspace)]) + desc    
					#    print desc_w, char_w, nspace
					#else:
					#    desc_w = len(desc[:last_space])
					#    nspace = int(np.floor((char_w - desc_w)/2))
					#    desc = desc[:last_space+1] + ''.join([' ' for _ in range(0, nspace)]) + desc[last_space+1:]
			#print desc, '+', x
			desc += x                                
			#print desc
                
		desc = '\n'.join([' ' + line.strip() for line in desc.split('\n')])
		#print raw_string(desc)
                

	box = ax_desc.get_window_extent()
	box_w, box_h = box.width, box.height

	desc_fmt = format_tex(desc)
	if text_obj is None:
		#print "START:", desc_fmt
		text_obj = ax_desc.text(
			0.5, 0.5, desc_fmt,
			va='top', ha='center',
			linespacing=1.25,fontsize=6, 
			bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1')
		)
	else:
		text_obj.set_text(desc_fmt)

	tbox = text_obj.get_window_extent(renderer=fig.canvas.get_renderer())
	tbox_w, tbox_h = tbox.width, tbox.height

	if char_w and char_w < 5:
		print "WARNING:", "char_w="+str(char_w)
		print desc
		print
		return False

	if tbox_w > box_w - 15:
		return add_description(fig, ax_desc, desc, text_obj=text_obj, char_w=char_w-1)
	#print "FINAL:", desc_fmt
	return True



###############################################################################
### file functions
###############################################################################
def generate_weblogo(inpath, target, outfile):
	""" Returns .weblogo.png image of sequence given by outfile """
	#print inpath
	#print target
	try:
		#return glob('/'.join([inpath, target, '*v2.png']))[0]
		return glob('/'.join([inpath, target, 'NATIVE-%s*.png' % target]))[0]
	except:
		return None


def save_figure(fig, fname, index=None, format='pdf'):
    ext = '.{fmt}'.format(fmt=format)
    fname = fname.replace(fname[-4:], ext)
    if index is not None:
        fname = fname.replace(
            ext, '.panel-{idx}{ext}'.format(idx=index, ext=ext)
        )
    fig.savefig(fname, format=format, dpi=1200)
    return fname

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
	order_file = options.order
	desc_file = options.desc
	debug = options.debug
	order = None
	if order_file:
		with open(order_file, 'r') as f:
			order = [l.strip() for l in f.readlines()[1:]]
	descriptions = parse_descriptions(desc_file)
            
	# check options
	inpaths = [abspath(x) for x in inpaths if exists(x) and isdir(x)]
	#assert(len(inpaths) == 1)

	#inpath = inpaths[0]
	inpath = None
       
	if targets[0] != '*':
		targets = targets
	elif target_files is not None:
		targets = get_target_names( target_files )
	else:
		targets = get_target_names( target_files, inpaths )

	if order:
		targets = order

	nplots = 15
	for target in targets:
		print "Target: "+target

                
        nplots = len(targets)
        print targets
        #if nplots >= 3:
        #        nplots = 3

        if pdfname is None:
			if inpath:
				pdfname = basename(inpath)+'_weblogos.pdf'
			else:
				pdfname = 'default_weblogos.pdf'
                

	# setup pdf and figure, return handles
	( pp, fullpdfname ) = setup_pdf_page( inpaths, targets, pdfname )
	#( fig, nplots, nrows, ncols ) = setup_figure( nplots, ncols = 3, landscape=False )
	( fig, nplots, nrows, ncols ) = setup_figure( nplots, landscape=False )
	fig.set_size_inches( 6.5, 1.5 * nrows)
	ncols = 3
	nrows = 3*int(np.ceil(nplots / ncols))+3
	
	fig.set_size_inches( 6.5, 1.5 * nrows)
        
	gridprops = dict(top=.95, bottom=.05,
                         left=.05, right=.95,
                         wspace=.2, hspace=.5)
	gs0 = gridspec.GridSpec(nrows, ncols)
	gs0.update(**gridprops)

	#ncols = 3
	#ncols *= 2
	nrows = int(np.ceil(nplots / ncols))
	nplots *= 2

        
	print ncols
	print nrows
	# iterate over targets
	count = 1
	panel_idx = 1
        
         
	for plot_idx, target in enumerate(targets, start=0):
                
		try:
			inpath = dirname(glob('*/'+target)[0])
		except:
			inpath = './'
                
		if target == 'BREAK' or count > nplots or count == 1 or plot_idx > len(targets):
			if count != 1:
				count = 1
				# get date printed to figure
				#plt.figtext(0.95,0.98,
				#            basename(inpath),
				#            horizontalalignment='right',
				#            fontsize='large')
				#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
				#lt.tight_layout()
				#fig.subplots_adjust(hspace=1.0)
				#fig.subplots_adjust(top=.95, bottom=.05,
				#                    left=.05, right=.95,
				#                    wspace=.3, hspace=1.5)

				#pp.savefig(DPI=1200)
				#plt.savefig(fullpdfname, DPI=1200)
				save_figure(fig, fullpdfname, index=panel_idx, format='svg')
				panel_idx += 1

			#( fig, nplots, nrows, ncols ) = setup_figure( nplots, ncols=3, landscape=False )
			( fig, nplots, nrows, ncols ) = setup_figure( nplots, landscape=False )
			#nplots = 12
			ncols = 3
			nrows = int(np.ceil(nplots / ncols))
			#ncols = 3
			print "nrows: ", nrows
			print "ncols: ", ncols
			gs0 = gridspec.GridSpec(nrows , ncols)
			#gs0.update(**gridprops)

			ncols *= 2
			nrows = int(np.ceil(nplots / ncols))
			#nplots *= 2

                        
		weblogo = generate_weblogo(inpath, target, silentfile)
		#print weblogo
                
		if weblogo is None or weblogo is False:
			print "target", target, "gave a bum image"
			continue
         
		gsX = gridspec.GridSpecFromSubplotSpec(
                        3, 1, 
                        height_ratios=[.25,1.0,.35],
						subplot_spec=gs0[count],#/2],
                        wspace=0.0, hspace=0.0
		) 
		#print gsX
                

		#############################################
		### Target Name
		#############################################
		#ax.set_title( target, x=1.0,y=1.25, fontsize=10, weight='bold' )
		ax_title = plt.subplot(gsX[0, :])
		ax_title.text(0.5, 0.5, format_tex(target, bold=True),
			verticalalignment='center', 
			horizontalalignment='center',
			fontsize=9, weight='bold',
			#bbox=dict(
			#         facecolor='none',
			#         edgecolor='none',
			#         boxstyle='round,pad=1'
			#)
		)
		plt.axis('off')         
		if debug:
			plt.axis('on')
                        
                
                
		#############################################
		### SS Graphic
		#############################################
		#ax_ss = plt.subplot(gsX[1, -1])

		# dummmy file, for testing
		#ss_filename = "/".join([os.path.dirname(__source_dir__), '/inputs/images/placeholder-secstruct.png'])
		#print
		#print "[warning] using place-holder secondary structure image... "
		#print "          ... need to implement automated generation/usage"
		#print "          using file: {}".format(ss_filename)
		#print 

		#ss_image = Image.open(ss_filename)
		#ss_image = ss_image.crop((0, 0, ss_image.size[0], ss_image.size[1]))
		#ss_image = ss_image.resize((ss_image.size[0], ss_image.size[1]))
		#img_arr = np.asarray(ss_image)

		#imgplt = plt.imshow(img_arr)#, aspect='equal')

		#plt.axis('off')               
		#if debug:
		#        plt.axis('on')
		#ax_ss.set_aspect('equal')
                
		#count += 1 
                

		#############################################
		### 3D Model
		#############################################
		# plot image on subplot
		ax = plt.subplot(gsX[1, 0])#, sharex=ax_ss, sharey=ax_ss) 
                
		image = Image.open(weblogo)
		image = image.crop((0, 0, image.size[0], image.size[1]))
		#r = float(ss_image.size[0])/float(image.size[0])
		#img_w = ss_image.size[0]
		#img_h = int(np.ceil(image.size[1]*r))
		#image = image.resize((img_w, img_h))
		img_arr = np.asarray(image)

		imgplt = plt.imshow(img_arr)#, aspect='equal')

		plt.axis('off')               
		if debug:  plt.axis('on')
		ax.set_aspect('equal')
                
		count += 1

		ax_desc = "NO DESCRIPTION FOUND"                

		#############################################
		### DESCRIPTION BOX
		#############################################
		if descriptions and target in descriptions:
			desc = descriptions[target]
			ax_desc = plt.subplot(gsX[-1, :])
			add_description(fig, ax_desc, desc)
			plt.axis('off')
			#if debug:
			#	plt.axis('on')
		else:
			desc = "Fill in later please"
			ax_desc = plt.subplot(gsX[-1, :])
			add_description(fig, ax_desc, desc)
			plt.axis('off')                                  
			#if debug:
			#	plt.axis('on')
                
		#for ax_ in [ax, ax_ss, ax_title]:#, ax_desc]:
		for ax_ in [ax, ax_title]:#, ax_desc]:
			ax_.set_xticklabels([])
			ax_.set_yticklabels([])
			for ticklabel in ax_.yaxis.get_ticklabels()+ax_.xaxis.get_ticklabels():
				ticklabel.set_fontsize(5)

	# finalize (adjust spacing, print date)
	#finalize_figure( fig, nplots, nrows, ncols )

	# save as pdf and close
	#plt.tight_layout(pad=0.4, w_pad=1.0, h_pad=1.0)
	#plt.tight_layout() #pad=4.0) #, h_pad=2.0, w_pad=0.4)
	#fig.subplots_adjust(hspace=1.0)
	#plt.subplots_adjust(top=.95, bottom=.1,
	#                    left=.08, right=.95,
	#                    wspace=.3, hspace=1.0)
	#plt.subplots_adjust(top=.95, bottom=.5,
	#                    left=.05, right=.95,
	#                    wspace=.3, hspace=1.0)
	#fig.subplots_adjust(top=.95, bottom=.05,
	#                    left=.05, right=.95,
	#                     wspace=.3, hspace=1.5)

	svgfilename = save_figure(fig, fullpdfname, index=panel_idx, format='svg')
	#svgfilename = save_figure(fig, fullpdfname, index=panel_idx, format='eps')
	plt.savefig(fullpdfname, DPI=1200) #pp.savefig(DPI=1200)
	#pp.close()

	# open pdf
	out, err = subprocess.Popen(['uname'], stdout=subprocess.PIPE).communicate()
	if 'Darwin' in out:
		subprocess.call(['open',svgfilename])
	if 'Linux' in out:
		subprocess.call(['xdg-open',svgfilename])

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
		'-order',
	    help='Ordered list of targets.',
	    default=None
	)
	
	parser.add_argument(
		'-o','--pdfname',
	    help='File name to save as pdf.',
	    default=None
	)
	parser.add_argument(
		'-desc',
	    help='File with descriptions of each target.',
	    default=None
	)
	parser.add_argument(
		'--debug',
	    help='File with descriptions of each target.',
	    action='store_true'
	)
	return parser


###############################################################################
### main
###############################################################################
if __name__=='__main__':

	sys.exit(make_weblogos(sys.argv[1:]))

