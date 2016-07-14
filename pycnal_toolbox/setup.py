#!/usr/bin/env python

"""
pycnal_toolbox is a suite of tools for working with ROMS.

Requires:
    pycnal (https://github.com/ESMG/PyCNAL)

Contains:
    many things...

"""

doclines = __doc__.split("\n")

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None,parent_package,top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True)
#                       quiet=True)
    config.add_subpackage('pycnal_toolbox')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "pycnal_toolbox",
          version = '0.1',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          author = "ESMG",
          url = 'https://github.com/ESMG/PyCNAL',
          license = 'BSD',
          platforms = ["any"],
          configuration=configuration,
          )
