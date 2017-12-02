import numpy as np
import pdb
import matplotlib.pyplot as plt
#import matplotlib.axes#.Axes# as axes

def main(filename):

    data = np.loadtxt(filename, delimiter=':') # h,m,s x 10000

    diffs = data[1:, :] - data[:-1,:] #time differences

    print "Mean", np.mean(diffs[:,2])
    print "Median", np.median(diffs[:,2])

    junk = plt.hist((diffs[:,2])**(-1), bins=200, log='T')
    #matplotlib.axes.Axes.set_ylim(0.5, 1000)
    plt.title('Distribution of SAPHIRA Frame Rates')
    plt.xlabel('Frame Rate')
    plt.ylabel('# of Occurences')
    plt.show()

    plt.plot(data[1:,2], diffs)
    plt.ylabel('Time difference between frames')
    plt.xlabel('Seconds')
    plt.show()
