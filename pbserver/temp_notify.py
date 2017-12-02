#!/usr/bin/env python

#Emails user if temp goes out of range. Needs to be restarted when SCExAO is 
# restarted.

import serial
import getpass
import datetime
import time
import smtplib
 
emailsent=0 #Has an email been sent yet?
testsent=0 #Has a test email been sent yet?

saph = serial.Serial('/dev/serial/by-id/usb-Silicon_Labs_Model_336_Temperature_Controller_336AATF-if00-port0', baudrate=57600, parity='O', bytesize=7, timeout=0.5)

server = smtplib.SMTP('smtp.gmail.com', 587)
server.starttls()
server.login('scexaonotifier@gmail.com', getpass.getpass())

while 1: #forever

    junk=saph.write("KRDG? A\r\n")
    A=saph.readlines()
    junk=saph.write("KRDG? B\r\n")
    B=saph.readlines()
    junk=saph.write("KRDG? C\r\n")
    C=saph.readlines()
    junk=saph.write("KRDG? D\r\n")
    D=saph.readlines()

    print "At " + str(datetime.datetime.now()) + " the GLS camera temperatures are:"
    print " A: " + str(A[0]) + " B: " + str(B[0]) + " C: " + str(C[0]) + " D: " + str(D[0])
    print

    #send a test email to verify code is working
    if testsent==0:
        msg = "The latest SAPHIRA temperatures are\n"+ \
              " A: " + str(A[0]) + " B: " + str(B[0]) + " C: " + str(C[0]) + " D: " + str(D[0])+'\n'+\
              "Congratulations, your code works.\n\n"\
              "Love,\n"\
              "Sean"
        server.sendmail("scexaonotifier@gmail.com", "geekyrocketguy@gmail.com", 
                        'Subject: Python Temperature Test Email\n'+msg)
        testsent=1


    #if temperature is out of range and user hasn't been notified yet
    if (float(A[0]) > 88) & (emailsent==0): 
        #send email
        msg = "The SAPHIRA temperatures are\n"+ \
              " A: " + str(A[0]) + " B: " + str(B[0]) + " C: " + str(C[0]) + " D: " + str(D[0])+'\n'+\
              "Thought you might want to know. Sorry about this news.\n\n"\
              "Love,\n"\
              "Sean"
        server.sendmail("scexaonotifier@gmail.com", "geekyrocketguy@gmail.com", 
                        'Subject: SAPHIRA Temp Warning\n'+msg)
        emailsent=1
        print "Email notification sent."

    time.sleep(300)

server.quit()
saph.close() #close serial connection
