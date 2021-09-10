import os.path as osp
from waflib.Configure import conf

def options(opt):
    opt = opt.add_option_group('VLNEval Options')

    opt.add_option(
        '--with-vlneval',
        type = 'string',
        help = "give VLNEval installation location",
    )

    opt.add_option(
        '--with-vlneval-include',
        default = '',
        help    = "give VLNEval include installation location",
        type    = 'string',
    )

    opt.add_option(
        '--with-vlneval-lib',
        default = '',
        help    = "give VLNEval lib installation location",
        type    = 'string',
    )

@conf
def check_vlneval(ctx, mandatory = True):
    instdir = ctx.options.with_vlneval

    if (instdir is None) or (instdir.lower() in ['yes', 'true', 'on']):
        ctx.start_msg('Checking for VLNEval in PKG_CONFIG_PATH')
        ctx.check_cfg(
            package      = 'vlneval',
            uselib_store = 'VLNEVAL',
            args         = '--cflags --libs',
            mandatory    = mandatory
        )

    elif instdir.lower() in ['no', 'off', 'false']:
        return

    else:
        ctx.start_msg('Checking for VLNEval in %s' % instdir)

        if ctx.options.with_vlneval_include:
            ctx.env.INCLUDES_VLNEVAL = [ ctx.options.with_vlneval_include ]
        else:
            ctx.env.INCLUDES_VLNEVAL = [ osp.join(instdir, 'include') ]

        if ctx.options.with_vlneval_lib:
            ctx.env.LIBPATH_VLNEVAL = [ ctx.options.with_vlneval_lib ]

        ctx.env.LIB_VLNEVAL = [ 'vlneval' ]

    ctx.check(
        header_name  = "vlneval/struct/VarDict.h",
        use          = 'VLNEVAL',
        uselib_store = 'VLNEVAL',
        define_name  = 'HAVE_VLNEVAL',
        mandatory    = mandatory
    )

    if len(ctx.env.INCLUDES_VLNEVAL):
        ctx.end_msg(ctx.env.INCLUDES_VLNEVAL[0])
    else:
        ctx.end_msg('VLNEVAL not found')

def configure(cfg):
    cfg.check_vlneval(mandatory = False)

    if cfg.env["HAVE_VLNEVAL"]:
        cfg.check_boost(lib = 'program_options')

