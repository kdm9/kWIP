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

from __future__ import print_function, division

from io import BytesIO

import pytest
from tempdir import run_in_tempdir

from kwipy.fastx import FastxReader, write_fastx


FASTA_SIMPLE = b"""\
>1 description
ACGTACGT
>2 description
ACGTACGT
"""
FASTA_NOEOL = b"""\
>1 description
ACGTACGT
>2 description
ACGTACGT"""
FASTA_MULTILINE = b"""\
>1 description
ACGT
ACGT
>2 description
ACGT
ACGT
"""

FASTA_BAD_HDR = b"""\
>1
"""


FASTA_PASS = [
    FASTA_SIMPLE,
    FASTA_NOEOL,
    FASTA_MULTILINE
]

FASTA_FAIL = [
]

FASTQ_SIMPLE = b"""\
@1 description
ACGTACGT
+
IIIIIIII
@2 description
ACGTACGT
+
IIIIIIII
"""

FASTQ_MULTILINE = b"""\
@1 description
ACGT
ACGT
+
IIII
IIII
@2 description
ACGT
ACGT
+
IIII
IIII
"""

FASTQ_PASS = [
    FASTQ_SIMPLE,
    FASTQ_MULTILINE
]

def test_read_fasta_good():
    for rec in FASTA_PASS:
        for i, seq in enumerate(FastxReader(BytesIO(rec), split_header=True)):
            assert seq.name == bytes(str(i + 1), 'ascii')
            assert seq.description == b'description'
            assert seq.sequence == b'ACGTACGT'
            assert seq.quality == b''
            assert seq.is_fastq == False

        for i, seq in enumerate(FastxReader(BytesIO(rec))):
            assert seq.name == bytes('{} description'.format(i + 1), 'ascii')
            assert seq.description == b''
            assert seq.sequence == b'ACGTACGT'
            assert seq.quality == b''
            assert seq.is_fastq == False

        # Wrong format
        with pytest.raises(ValueError):
            for seq in FastxReader(BytesIO(rec), format='fastq'):
                pass


def test_read_fasta_bad():
    for rec in FASTA_FAIL:
        with pytest.raises(ValueError):
            for seq in FastxReader(BytesIO(rec)):
                pass


def test_read_fastq_good():
    for rec in FASTQ_PASS:
        for i, seq in enumerate(FastxReader(BytesIO(rec), split_header=True)):
            assert seq.name == bytes(str(i + 1), 'ascii')
            assert seq.description == b'description'
            assert seq.sequence == b'ACGTACGT'
            assert seq.quality == b'IIIIIIII'
            assert seq.is_fastq == True

        for i, seq in enumerate(FastxReader(BytesIO(rec))):
            assert seq.name == bytes('{} description'.format(i + 1), 'ascii')
            assert seq.description == b''
            assert seq.sequence == b'ACGTACGT'
            assert seq.quality == b'IIIIIIII'
            assert seq.is_fastq == True
