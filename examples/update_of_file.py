#!/usr/bin/env python


file = open('/tmp/workfile', 'r')
import time
while 1:
    where = file.tell()
    line = file.readline()
    if not line:
        time.sleep(1)
        file.seek(where)
    else:
        y.append(float(line))
