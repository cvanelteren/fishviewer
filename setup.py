# from distutils.core import setup, Extension
from setuptools import setup
from setuptools.extension import Extension
# MODULE INFORMATION
__module__              = "fishviewer"
__author__              = "Casper van Elteren"
__email__               = "caspervanelteren@gmail.com"
__description__         = "Fast data viewer for zebrafish data"
__license__             = "MIT"
__python__requires__    = ">=3.6"


with open('requirements.txt', 'r') as f:
    requirements =  f.readlines()
requirements = [i.strip() for i in requirements]
# compile
with open('requirements.txt', 'r') as install_reqs:
    setup(\
        name            = __name__,\
        author          = __author__,\
        author_email    = __email__,\
        python_requires = __python__requires__,\
        zip_safe        = False,\
        install_requires= requirements,\

    # gdb_debug =True,
    )
