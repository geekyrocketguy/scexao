#!/usr/bin/env python

#Emails user if temp goes out of range. Is called with a chronjob.

#Google will occasionally get grumpy and block the login attempt (maybe
# if the IP address changes?). Either open a browser on SCExAO2 and log
# in to scexaonotifier@gmail.com, or click the "Yes that was me" on the
# "Someone has your password!" security emails that gmail sends out.
# Then everything will be happy again.

import serial
#import getpass
#import datetime
#import time
import smtplib
import os.path
 
temp_threshold = 87 #At what temp should the user be notified?
to_address = 'geekyrocketguy@gmail.com' #Who should the email be sent to?

#emailsent=0 #Has an email been sent yet?
#testsent=0 #Has a test email been sent yet?


saph = serial.Serial('/dev/serial/by-id/usb-Silicon_Labs_Model_336_Temperature_Controller_336AATF-if00-port0', baudrate=57600, parity='O', bytesize=7, timeout=0.5)

junk=saph.write("KRDG? A\r\n")
A=saph.readlines()
junk=saph.write("KRDG? B\r\n")
B=saph.readlines()
junk=saph.write("KRDG? C\r\n")
C=saph.readlines()
junk=saph.write("KRDG? D\r\n")
D=saph.readlines()
saph.close() #close serial connection

#UPDATE 4 11 2018
#B = [A[1]]
#C = [A[2]]
#D = [A[4]]
#A = [A[0]]

if not os.path.isfile('emailed.txt'): #if someone deleted the file, recreate it
    f=open('emailed.txt', 'w')
    f.write('0')
    f.close()
    print "emailed.txt was deleted by some goon, but it has been restored."

f=open('emailed.txt', 'r')
emailsent=int(f.read()) #has the user been emailed recently?
f.close()

if emailsent: #if the user has already been emailed:
    if (float(A[0]) < (temp_threshold-1) ): #if temperature is nearly back in range
        f=open('emailed.txt', 'w') #reset the status so the user can be emailed again
        f.write('0')
        f.close()

        #send email saying the detector is cooled again and ready to use
        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.starttls()
        pw = np.loadtxt('notthepassword.txt', dtype='str')
        server.login('scexaonotifier@gmail.com', pw)
        msg = "The SAPHIRA temperatures are\n"+ \
              " A: " + str(A[0]) + " B: " + str(B[0]) + " C: " + str(C[0]) + " D: " + str(D[0])+'\n'+\
              "Thought you might want to know. Congratulations on this news.\n\n"\
              "Love,\n"\
              "Sean"
        server.sendmail("scexaonotifier@gmail.com", to_address, 
                        'Subject: SAPHIRA Temp Back to Normal\n'+msg)
        server.quit()
 
    else: #if the temperature is still out of range, do nothing
        print "SAPHIRA temperature is out of range, but user has already been emailed about this."

else: #if the user has not been emailed
    if (float(A[0]) > temp_threshold): #if temperature is out of range, email the user!

        #send email
        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.starttls()
        server.login('scexaonotifier@gmail.com', '')
        msg = "The SAPHIRA temperatures are\n"+ \
              " A: " + str(A[0]) + " B: " + str(B[0]) + " C: " + str(C[0]) + " D: " + str(D[0])+'\n'+\
              "Thought you might want to know. Sorry about this news.\n\n"\
              "Love,\n"\
              "Sean"
        server.sendmail("scexaonotifier@gmail.com", to_address, 
                        'Subject: SAPHIRA Temp Warning\n'+msg)
        server.quit()
        
        #don't email the user again
        f=open('emailed.txt', 'w')
        f.write('1')
        f.close()

    else:#, if the temperature is in range, do nothing
        print "SAPHIRA temperature is in range, everything is happy."

