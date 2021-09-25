#!/usr/bin/env python


TOP = '.'
APPNAME = 'Junk'

def options(opt):
    opt.load("wcb")

def configure(cfg):
    cfg.load("wcb")

    if cfg.env["HAVE_VLNEVAL"]:
        cfg.check_boost(lib = 'program_options')

    # boost 1.59 uses auto_ptr and GCC 5 deprecates it vociferously.
    cfg.env.CXXFLAGS += ['-Wno-deprecated-declarations']
    cfg.env.CXXFLAGS += ['-Wall', '-Wno-unused-local-typedefs', '-Wno-unused-function']
    # cfg.env.CXXFLAGS += ['-Wpedantic', '-Werror']


def build(bld):
    bld.load('wcb')

    use = [ 'ROOTSYS' ]

    if bld.env["HAVE_VLNEVAL"]:
        use += [ 'VLNEVAL', 'BOOST' ]

    bld.smplpkg('WCPuBooNE_BDT_APP', use = use)

