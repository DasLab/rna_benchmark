import os, errno
from os import chdir,getcwd
from os.path import basename,dirname,exists
import glob
from sys import argv
import subprocess

CWD = getcwd()
chdir( dirname( argv[0] ) )

if not exists('bin'):
    subprocess.call(['mkdir','bin'])

for f in glob.glob('./bin/*'):
    os.remove(f)

for dirpath, dirnames, filenames in os.walk('./scripts/'):
    if len( basename(dirpath) ) == 0 or basename(dirpath)[0] == '.' or basename(dirpath) == 'bin':
        continue
    for f in filenames:
        filename = os.path.join(dirpath, f)
        if filename[-3:] == '.py':
            #print filename
            file1 = '../' + filename
            file2 = './bin/' + f
            try:
                os.symlink(file1, file2)
            except OSError, e:
                if e.errno == errno.EEXIST:
                    os.remove(file2)
                    os.symlink(file1, file2)

chdir( CWD )
