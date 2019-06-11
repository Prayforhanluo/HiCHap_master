# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 10:06:52 2018

@author: han-luo
"""

import sys, HiCHap, os, glob, subprocess
from setuptools import  setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if (sys.version_info.major != 2) or (sys.version_info.minor < 6):
    print 'PYTHON VERSION MUST BE 2.6 or 2.7. YOU ARE CURRENTLY USING PYTHON ' + sys.version
    sys.exit(2)

# Guarantee Unix Format
for src in glob.glob('scripts/*'):
    text = open(src, 'r').read().replace('\r\n', '\n')
    open(src, 'w').write(text)
    
setup(name = 'HiCHap',
      version = HiCHap.__version__,
      author = HiCHap.__author__,
      url = 'https://github.com/Prayforhanluo/HiCHap_master',
      author_email = 'hluo_lc@outlook.com',
      description = 'A Library of Hi-C data processing, bias correction and structural analysis for phased haplotype',
      long_description = read('README.rst'),
      long_description_content_type = 'text/x-rst',
      scripts = glob.glob('scripts/*'),
      packages = find_packages(),
      classifiers = [
        'Programming Language :: Python',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)']
      )
      

print "Fix the pkg_resources.ResolutionError::"

execute = subprocess.Popen('which hichap', shell = True, stdout=subprocess.PIPE)

path = execute.stdout.read().strip()
print "Executable file path  : %s " % path
print "Copying build/scripts-2.7/hichap -> %s" %path

subprocess.call('cp ./build/scripts-2.7/hichap %s' % path, shell = True)

print "Install Done!"
