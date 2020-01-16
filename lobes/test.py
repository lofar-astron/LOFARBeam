#!/usr/bin/env python3

import numpy as np

import lobes

a = np.linspace(0, 0.5, 10)

lobes.scale(a, 2.0)

print(lobes.legendre(1,-11,a))

