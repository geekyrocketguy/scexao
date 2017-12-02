#!/usr/bin/env python

#Takes as input a string. Finds the thread ID in it.

import numpy as np
import pdb
import sys
   
args=sys.argv

if len(args) != 2: 
    print
    print "Usage: python find_tid.py string"
    print "Returns the thread ID."
    print
    #return
else:
    mystr=str(args[1])[42:] #leaves only PID and TID
    print mystr[mystr.find(':')+1:]
