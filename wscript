#! /usr/bin/env python
import sys
import os
import sferes
sys.path.insert(0, sys.path[0]+'/waf_tools')
print sys.path[0]

#import limbo

from waflib.Configure import conf

def build(bld):
    


    bld.program(features = 'cxx',
                source = 'test_eval.cpp',
                includes = '. ../../ ../../src /git/sferes2/',
                uselib = 'TBB BOOST EIGEN PTHREAD MPI LIMBO',
                use = 'limbo',
                target = 'test_eval')

def configure(conf):
    pass

def options(opt):
    pass

