from setuptools import setup, find_packages
from distutils.core import Extension
from distutils import sysconfig
import numpy as np
import platform
import os

#Hack to get rid of the excess compiler flags distutils adds by default...
#Right now I just strip the -g flag 
#http://stackoverflow.com/questions/13143294/distutils-how-to-disable-including-debug-symbols-when-building-an-extension
#Remove the code below if it causes problems and you don't care about using the debug flag
if platform.system() != 'Windows':  # When compilinig con visual no -g is added to params
    cflags = sysconfig.get_config_var('CFLAGS')
    opt = sysconfig.get_config_var('OPT')
    sysconfig._config_vars['CFLAGS'] = cflags.replace(' -g ', ' ')
    sysconfig._config_vars['OPT'] = opt.replace(' -g ', ' ')

if platform.system() == 'Linux':  # In macos there seems not to be -g in LDSHARED
    ldshared = sysconfig.get_config_var('LDSHARED')
    sysconfig._config_vars['LDSHARED'] = ldshared.replace(' -g ', ' ')


# An extension configuration has the following format:
# module_name' : { 'extension kwarg' : argument }
# This makes adding C from deep in the package trivial.  No more nested setup.py.

#I can't get the wildcard for the source .c files to work for some reason. 
source_path = os.path.join(os.path.dirname(__file__), 'surfwave_inv', 'sk_disp_crv/')
extension_config = {'surfwave_inv.sk_disp_crv._sk_disp_crv' :
                          { 'sources' : [source_path +'get_disp_crv.c', source_path +'disp_fun.c'],
                            'extra_compile_args' :  ["-O2","-ffast-math"], #FOR SOME REASON I RUN INTO PROBLEMS WITH -03 SOMETIMES
                            'include_dirs' : [np.get_include(), os.path.join(os.path.dirname(__file__), 'surfwave_inv','sk_disp_crv')]
                          },
                   }

extensions = [Extension(key, **value) for key, value in extension_config.iteritems()]

setup( # Update these.  Replace 'template' with the name of your extension.
    version = '0.0', 
    name = 'surfwave_inv_pkg',
    description = 'A wrapped C code for evaluating dispersion curves using Schwab and Knopoff 1972 applied to an exploration geophysics setting.',
    ext_modules = extensions,
    
    # Do not change any of these unless you know what you are doing.
    install_requires = ['setuptools'],
    packages=find_packages(),
    namespace_packages = ['surfwave_inv'],
    zip_safe = False
    )
