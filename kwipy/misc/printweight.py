#!/usr/bin/env python3
import bcolz
import numpy as np

import sys

print(bcolz.open(sys.argv[1])[:])
