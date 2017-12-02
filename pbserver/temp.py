#!/usr/bin/env python
import serial

saph = serial.Serial('/dev/serial/by-id/usb-Silicon_Labs_Model_336_Temperature_Controller_336AATF-if00-port0', baudrate=57600, parity='O', bytesize=7, timeout=0.5)

junk=saph.write("KRDG? A\r\n")
A=saph.readlines()
junk=saph.write("KRDG? B\r\n")
B=saph.readlines()
junk=saph.write("KRDG? C\r\n")
C=saph.readlines()
junk=saph.write("KRDG? D\r\n")
D=saph.readlines()

print "The GLS camera temperatures are:"
print " A: " + str(A[0]) + \
    " B: " + str(B[0]) + \
    " C: " + str(C[0]) + \
    " D: " + str(D[0])

saph.close()

#junk=saph.write("SETP 1,85.000\r\n") #set A to 85K (will ramp down)
#junk=saph.write("SETP? A\r\n") #check setpoint of A
