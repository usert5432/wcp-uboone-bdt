import os.path as osp
from waflib.Configure import conf

def options(opt):
    opt = opt.add_option_group('VLNetEval Options')

    opt.add_option(
        '--with-vlneval',
        type = 'string',
        help = "give VLNetEval installation location",
    )

    opt.add_option(
        '--with-vlneval-include',
        default = '',
        help    = "give VLNetEval include installation location",
        type    = 'string',
    )

    opt.add_option(
        '--with-vlneval-lib',
        default = '',
        help    = "give VLNetEval lib installation location",
        type    = 'string',
    )

@conf
def check_vlneval(ctx, mandatory = True):
    ctx.end_msg('fdfdsafdsfdasfdasfdas')
    instdir = ctx.options.with_vlneval

    if (instdir is None) or (instdir.lower() in ['yes', 'true', 'on']):
        ctx.start_msg('Checking for VLNetEval in PKG_CONFIG_PATH')
        ctx.check_cfg(
            package      = 'vlneval',
            uselib_store = 'VLNETEVAL',
            args         = '--cflags --libs',
            mandatory    = mandatory
        )

    elif instdir.lower() in ['no', 'off', 'false']:
        return

    else:
        ctx.start_msg('Checking for VLNEval in %s' % instdir)

        if ctx.options.with_vlneval_include:
            ctx.env.INCLUDES_VLNETEVAL = [ ctx.options.with_vlneval_include ]
        else:
            ctx.env.INCLUDES_VLNETEVAL = [ osp.join(instdir, 'include') ]

        if ctx.options.with_vlneval_lib:
            ctx.env.LIBPATH_VLNETEVAL = [ ctx.options.with_vlneval_lib ]

        ctx.env.LIB_VLNETEVAL = [ 'vlneval' ]

    ctx.check(
        header_name = "vlneval/struct/VarDict.h",
        use         = 'VLNETEVAL',
        mandatory   = mandatory
    )

    if len(ctx.env.INCLUDES_VLNETEVAL):
        ctx.end_msg(ctx.env.INCLUDES_VLNETEVAL[0])
    else:
        ctx.end_msg('VLNETEVAL not found')

def configure(cfg):
    cfg.check_vlneval()
    cfg.check_boost(lib='program_options')

