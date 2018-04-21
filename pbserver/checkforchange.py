#!/usr/bin/env python

#Checks site to see if it has changed.
#syntax: python checkforchange.py 20180709 test
#date is required, test is optional

#Google will occasionally get grumpy and block the login attempt (maybe
# if the IP address changes?). Either open a browser on SCExAO2 and log
# in to scexaonotifier@gmail.com, or click the "Yes that was me" on the
# "Someone has your password!" security emails that gmail sends out.
# Then everything will be happy again.

#import serial
#import getpass
#import datetime
#import time
import smtplib
import os.path
import urllib
#import numpy as np
import sys
#import pdb

args=sys.argv
date = args[1]
 
to_address = 'geekyrocketguy@gmail.com' #Who should the email be sent to?
url='https://camping.ehawaii.gov/camping/all,sites,0,25,1,1692,UNDESIGNATED,,,'+date+',5,,,,1,.html'

page = urllib.urlopen(url)
pagecontents = page.read()
pagecontents = pagecontents[pagecontents.find('Campground Type') : pagecontents.find('Milolii')]

#clean HTML
while pagecontents.find('<') != -1: #remove html tags
    pagecontents = pagecontents.replace(pagecontents[pagecontents.find('<') : pagecontents.find('>')+1], '')
pagecontents = pagecontents.replace('&nbsp;', '')
pagecontents = pagecontents.replace(' ', '') #remove excess spaces
while pagecontents.find('\n\n') != -1: #remove excess line spacing
    pagecontents = pagecontents.replace('\n\n', '\n')

#check if file exists
if not os.path.isfile('avail_'+date+'.txt'): #if someone deleted the file, recreate it
    f=open('avail_'+date+'.txt', 'w')
    f.write(pagecontents)
    f.close()
    print 'avail_'+date+'.txt'+" was deleted by some goon, but it has been restored."

f=open('avail_'+date+'.txt', 'r')
oldcontents=f.read() #has the user been emailed recently?
f.close()

if (pagecontents != oldcontents) or ('test' in args): #has something changed? Then email user.
        #send email saying the detector is cooled again and ready to use
        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.starttls()
        pw = np.loadtxt('notthepassword.txt', dtype='str')
        server.login('scexaonotifier@gmail.com', pw)

        if 'test' in args:
            msg = 'THIS IS A TEST\n\n'
            mysubject = 'Code is Working'
        else:
            msg = ''
            mysubject = 'Change in Permit Availability'

        msg += "The campground availability is as follows:\n" + \
               pagecontents + \
               "Thought you might want to know.\n\n"\
               "Love,\n"\
               "Sean"
        
        server.sendmail("scexaonotifier@gmail.com", to_address, 
                        'Subject: '+mysubject+'\n'+msg)
        server.quit()

        #print new availability into text document
        f=open('avail_'+date+'.txt', 'w')
        f.write(pagecontents)
        f.close()
else:
    print "Nothing has changed, code is happy."
