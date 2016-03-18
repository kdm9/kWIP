#!/usr/bin/env python

# Copyright  2016  Kevin Murray <spam@kdmurray.id.au>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from setuptools import setup
from setuptools.command.test import test as TestCommand
EXT='c'
try:
    from Cython.Build import cythonize
    EXT='pyx'
except ImportError:
    def cythonize(x): return x

import versioneer


class NoseCommand(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        # Run nose ensuring that argv simulates running nosetests directly
        import nose
        nose.run_exit(argv=['nosetests'])


desc = """
kwip in python
"""

setup_requires = [
    'nose',
    'cython',
]

with open('requirements.txt') as fh:
    install_requires = [req.strip() for req in fh]

test_requires = [
    'coverage>=3.7',
]

command_classes=versioneer.get_cmdclass()
command_classes['test'] =  NoseCommand

setup(
    name="kwipy",
    packages=['kwipy', ],
    version=versioneer.get_version(),
    entry_points={
        'console_scripts': [
            'kwipy-count = kwipy.scripts:count_main',
            'kwipy-weight = kwipy.scripts:weight_main',
            'kwipy-kernel-mpi = kwipy.scripts:kernel_mpi_main',
            'kwipy-distmat = kwipy.scripts:distmat_main',
        ],
    },
    cmdclass=command_classes,
    install_requires=install_requires,
    tests_require=test_requires,
    setup_requires=setup_requires,
    description=desc,
    author="Kevin Murray",
    author_email="spam@kdmurray.id.au",
    url="https://github.com/kdmurray91/kwip",
    keywords=[
        "bioinformatics",
        "Next-gen Sequencing",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
)
