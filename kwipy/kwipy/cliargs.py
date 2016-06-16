# Copyright 2016 Kevin Murray <spam@kdmurray.id.au>
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

from __future__ import print_function, division, absolute_import

import argparse
from argparse import (
    ArgumentParser,
    RawDescriptionHelpFormatter as RawFormatter,
)
from textwrap import dedent


def count_args():
    desc = dedent("""\
    Counts k-mers into individual count vectors, parallelised using MPI.""")
    epilog = dedent('''\
    Will use about 6 * CVLEN bytes of RAM per file (or 2x with --no-cms).

    An optional pre-counting command for e.g. QC or SRA dumping can be given
    with `-p`. The pre-command can be a shell pipeline combining the effects of
    multiple programs. Interleaved or single ended reads must be printed on
    stdout by the command(s). The pre-command uses the find/xargs/GNU Parallel
    convention of using a pair of '{}' to mark where the filename should be
    placed. Examples of a pre-command include:

        --precmd 'fastq-dump --split-spot --stdout {}'
        --precmd 'zcat {} | trimit'
    ''')

    parser = ArgumentParser(description=desc, epilog=epilog,
                            formatter_class=RawFormatter)

    parser.add_argument(
        '-k', '--ksize', type=int, default=20,
        help='K-mer length')
    parser.add_argument(
        '-v', '--cvsize', type=float, default=5e8,
        help='Count vector length')
    parser.add_argument(
        '-p', '--precmd', required=False,
        help='Shell pipeline to run on input files before hashing')
    parser.add_argument(
        '--no-cms', action='store_false', dest='use_cms',
        help='Disable the CMS counter, use only a count vector')
    parser.add_argument('outfile', help='Output file/directory')
    parser.add_argument('readfiles', nargs='+', help='Read files')
