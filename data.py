#!/usr/bin/env python3

import json

with open('data.json', 'r') as f:
    samples = json.load(f)

print('#ifndef _DATA_H_')
print('#define _DATA_H_')
print('#define STATE_N 4')
print('#define DATA_N {}'.format(len(samples)))
print('#define NODE_N {}'.format(len(samples[0])))
print('int data[] = {')
for sample in samples:
    print(','.join(map(lambda b: str(int(b)), sample.values())) + ',')
print('};')
print('#endif')
