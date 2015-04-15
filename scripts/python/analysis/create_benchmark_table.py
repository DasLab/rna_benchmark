#!/usr/bin/python

################################################################################
### import modules
################################################################################
import os
import argparse
import subprocess as sp
import multiprocessing as mp
from os.path import exists, basename, dirname, isdir, isfile, realpath, abspath
from glob import glob
from create_benchmark_table_util import *


################################################################################
### main function
################################################################################
def get_target_row(target, (args)):
    
    ############################################################################
    ### change into target_dir and find all silent files
    ############################################################################
    working_dir = os.getcwd()
    os.chdir(target)    
    
    ############################################################################
    ### get table info
    ############################################################################
    table_row = TableRow()
    table_row.add_columns( get_title( target ) )
    table_row.add_columns( get_target_properties( target ) )
    table_row.add_columns( get_best_of_lowest_energy_cluster_centers( target, args.silent_file_name ) )
    table_row.add_columns( get_lowest_rmsd_model( args.silent_file_name ) )
    table_row.add_columns( get_lowest_energy_sampled( args.silent_file_name ) )

    ############################################################################
    ### change back to working directory
    ############################################################################
    os.chdir(working_dir)
    return table_row


################################################################################
### main script
################################################################################
if __name__=='__main__':


    ############################################################################
    ### parse arguments
    ############################################################################
    parser = argparse.ArgumentParser(
        description='Create a table of lowest energies/rmsds, \
                     and best of 5 cluster center rmsds from \
                     a list of runs (one table per run).'
    )

    parser.add_argument(
        'inpaths',
        nargs='+',
        help='List of paths to target runs.'
    )
    parser.add_argument(
        '-t','--targets',
        nargs='+',
        help='List of targets.',
        default=[]
    )
    parser.add_argument(
        '-s','--silent_file_name',
        help='Names of final silent files.',
        default=['region_FINAL.out','swm_rebuild.out']
    )
    parser.add_argument(
        '-f','--force',
        help='Rewrite tables even if they alread exist.',
        action='store_true'
    )
    parser.add_argument(
        '-j','--nproc',
        help='Number of jobs to run in parallel',
        default=(mp.cpu_count()-1)
    )

    args = parser.parse_args()
    inpaths = sorted( args.inpaths )
    user_targets = sorted( args.targets )
    force = args.force
    nproc = int(args.nproc)


    ############################################################################
    ### checks and initializations 
    ############################################################################
    assert( all( exists( inpath ) for inpath in inpaths ) )
    working_dir = os.getcwd()

    ############################################################################
    ### change into inpaths and get targets
    ############################################################################
    for idx, inpath in enumerate(inpaths, start=1):
        
        os.chdir( inpath )
        print '\n[%d/%d] Creating Table for Run: %s' % (idx,len(inpaths),inpath)

        ########################################################################
        ### find all targets in inpath
        ########################################################################
        found_targets = sorted(filter(isdir, glob('*')))
        if len(user_targets):
            found_targets = filter(lambda x: x in user_targets, found_targets)
        targets = filter(lambda x: x in found_targets, get_target_names())


        ########################################################################
        ### get table info for all targets found in inpath
        ########################################################################
        if nproc > 0:
            pool = mp.Pool(processes=nproc)

            out = [pool.apply_async(get_target_row,args=(t,args)) for t in targets]
            table_row_list = [o.get() for o in out]
            
            pool.close()
            pool.join()

        else:
            table_row_list = [get_target_row(t,args) for t in targets]

        ########################################################################
        ### create table for inpath
        ########################################################################
        table_name = inpath.upper() + '.tbl'
        if exists( table_name ) and not force:
            print "Table:", table_name, "already exists!!!"
        else:
            table = Table( table_name )
            table.add_row( column_names )
            table.add_row([' | '.join(subcolumn_names[c]) for c in column_names])
            for table_row in table_row_list:
                table.add_row( table_row.get_columns() )
            table.save()

        ########################################################################
        ### change back into working directory
        ########################################################################
        os.chdir( working_dir )

