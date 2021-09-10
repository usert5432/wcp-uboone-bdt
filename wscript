#!/usr/bin/env python


TOP = '.'
APPNAME = 'Junk'

def options(opt):
    opt.load("wcb")
    opt.load('boost')
    opt.load("find_vlneval")

def configure(cfg):
    cfg.load("wcb")
    cfg.load('boost')
    cfg.load("find_vlneval")

    cfg.write_config_header("config.h")

    # boost 1.59 uses auto_ptr and GCC 5 deprecates it vociferously.
    cfg.env.CXXFLAGS += ['-Wno-deprecated-declarations']
    cfg.env.CXXFLAGS += ['-Wall', '-Wno-unused-local-typedefs', '-Wno-unused-function']
    # cfg.env.CXXFLAGS += ['-Wpedantic', '-Werror']


def build(bld):
    bld.load('wcb')
    bld.load("find_vlneval")

    use = [ 'ROOTSYS' ]

    if bld.env["HAVE_VLNEVAL"]:
        use += [ 'VLNEVAL', 'BOOST' ]

    bld.smplpkg('WCPuBooNE_BDT_APP', use = use)

