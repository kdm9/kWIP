from __future__ import print_function, division, absolute_import
from kwipy.progress import ProgressLogger
from time import time, sleep
from random import random

bar = ProgressLogger(interval=10)
for i in range(100):
    sleep(0.1*random())
    bar.update(i + 1)
bar.done()

print("As a wrapper")
for i in ProgressLogger(range(104), 'reads', interval=10):
    sleep(0.02*random())
