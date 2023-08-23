#!/bin/env python3

import sys
import numpy as np

def polynomial(x, cs):
    r = 0
    for c in cs:
        r = r * x + c
    return r

point_count = 5000

coefficients = {
        'point': [0],
        'line1': [1, 0],
        'line2': [2, 0],
        'line3': [0.2, 0],
        'arc1': [-1, 0, 0],
        'arc2': [2, 1, 0],
        'arc3': [1, -2, 1],
        'es1': [1, 0, -1, 0],
        'es2': [1, 0, -2, 0],
        'es3': [-0.5, 0, 2, 0]
    }


for name, cs in coefficients.items():
    for q in range(1, 30):
        data = []
        for i in range(0, point_count):
            x = 0
            y = 0
            if len(cs) > 1:
                x = np.random.random() * 4 - 2
                y = polynomial(x, cs)

            x += np.random.normal(0, 0.015 * q)
            y += np.random.normal(0, 0.015 * q)

            data.append((x, y))
        with open("data/{}_{}.data".format(name, q), "w") as f:
            for x, y in data:
                f.write("{},{}\n".format(x, y))
