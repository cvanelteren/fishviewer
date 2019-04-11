# from distutils.core import setup, Extension
from setuptools import setup, find_packages
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
print(requirements)
setup(\
    version         = "1.0",\
    name            = __name__,\
    author          = __author__,\
    author_email    = __email__,\
    python_requires = __python__requires__,\
    packages        = find_packages(),\
)
