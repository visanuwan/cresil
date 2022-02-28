import os
import sys
import re
from subprocess import check_call
from setuptools import setup, find_packages
from setuptools.command.install import install


__pkg_name__ = 'cresil'

verstrline = open(os.path.join(__pkg_name__, '__init__.py'), 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}/__init__.py".'.format(__pkg_name__))

with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="{}".format(__pkg_name__),
    version=__version__,
    packages=find_packages(),
    include_package_data=True,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Visanu Wanchai, Piroon Jenjaroenpun, Intawat Nookaew',
    author_email='visanu.wanchai@gmail.com',
    url='https://gitlab.com/visanu/cresil',
    entry_points = {
        'console_scripts': [
            '{0} = {0}:main'.format(__pkg_name__)
        ]
    },
)