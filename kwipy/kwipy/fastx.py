from collections import namedtuple
from io import BufferedReader
import gzip
import bz2
import lzma

SeqRecord = namedtuple("SeqRecord", ['name', 'description', 'sequence',
                                     'quality', 'is_fastq'])

FASTA = 1
FASTQ = 2

ZIP_FMTS = [
    (b'\x1f\x8b', gzip.open),
    (b'BZ', bz2.open),
    (b'\xfd7zXZ\x00', lzma.open),
]


def decompress(stream):
    if not hasattr(stream, "peek"):
        stream = BufferedReader(stream)
    firstbytes = stream.peek(10)
    for magic, zwrapper in ZIP_FMTS:
        if firstbytes[:len(magic)] == magic:
            return zwrapper(stream, 'rb')
    return stream


class FastxReader(object):
    '''Iterator over a fasta/q records in some byte stream.'''

    def __init__(self, stream, format=None, split_header=False):
        self.stream = decompress(stream)
        if isinstance(format, str):
            if format.lower() == 'fasta':
                self.format = FASTA
            elif format.lower() == 'fastq':
                self.format = FASTQ
            else:
                raise ValueError("invalid format= arg")
        else:
            self.format = format
        self.split_header = split_header
        self.line = None

    def _check_format(self, line):
        if self.format is None:
            if line.startswith(b'>'):
                self.format = FASTA
            elif line.startswith(b'@'):
                self.format = FASTQ
            else:
                raise ValueError("Malformatted input stream, cannot parse "
                                 "FASTA/Q record")
        else:
            if ((self.format == FASTA and not line.startswith(b'>')) or
                    (self.format == FASTQ and not line.startswith(b'@'))):
                raise ValueError("Malformatted input stream, cannot parse "
                                 "FASTA/Q record")

    def _fasta(self):
        sequence = []
        for line in self.stream:
            if line.startswith(b'>'):
                self.line = line
                break
            sequence.append(line.strip())
        return b''.join(sequence)

    def _fastq(self):
        sequence = []
        quality = None
        for line in self.stream:
            if quality is None:
                if line.startswith(b'+'):
                    quality = []
                    continue
                else:
                    sequence.append(line.rstrip())
            elif len(quality) < len(sequence):
                quality.append(line.rstrip())
            else:
                self.line = line
                break
        sequence = b''.join(sequence)
        quality = b''.join(quality)
        if len(sequence) != len(quality):
            raise ValueError("Invalid fastq stream")
        return sequence, quality

    def __iter__(self):
        return self

    def __next__(self):
        name = description = sequence = quality = b''
        while True:
            # Get line
            if self.line is not None:
                line = self.line
                self.line = None
            else:
                line = self.stream.readline().rstrip()
                if not line:
                    raise StopIteration

            # Set format
            self._check_format(line)

            # Header
            line = line.rstrip()[1:]
            if self.split_header:
                name, _, description = line.partition(b' ')
            else:
                name = line

            # Get seq and qual
            if self.format == FASTA:
                sequence = self._fasta()
            elif self.format == FASTQ:
                sequence, quality = self._fastq()

            # make & yield record
            is_fastq = self.format == FASTQ
            return SeqRecord(name, description, sequence, quality, is_fastq)


def write_fastx(seqrec, wrap=80):
    pass
