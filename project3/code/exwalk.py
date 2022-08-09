# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:19:57 2013

@author: svenni
"""

from pylab import *
from scipy.ndimage import measurements
from walk import walk

# exwalk.py
# Generate spanning cluster (l-r spanning)

lx = 100
ly = 100
p = 0.586
ncount = 0
perc = []
while (len(perc)==0):
    ncount = ncount + 1
    if (ncount >1000):
        print("Couldn't make percolation cluster...")
        break

    z=rand(lx,ly)<p
    lw,num = measurements.label(z)
    perc_x = intersect1d(lw[0,:],lw[-1,:])
    perc = perc_x[where(perc_x > 0)]
    print(ncount)

if len(perc) > 0:
    labelList = arange(num + 1)
    area = measurements.sum(z, lw, index=labelList)
    areaImg = area[lw]
    maxArea = area.max()
    zz = (lw == perc[0])

    # zz now contains the spanning cluster
    figure()
    title("Spanning cluster")
    imshow(zz, interpolation='nearest', origin='upper') # Display spanning cluster
    # savefig("current.pdf")

    #% Run walk on this cluster
    l,r = walk(zz)
    figure()
    title("Left walker")
    imshow(l, interpolation='nearest', origin='upper')
    figure()
    title("Right walker")
    imshow(r, interpolation='nearest', origin='upper')
    zzz = l*r # Find points where both l and r are non-zero

    figure()
    title("Single connected points")
    imshow(zzz, interpolation='nearest', origin='upper')
    colorbar()

    M_SC = np.sum(zzz > 0)
    print(M_SC)


    # print(M_SC)
    #print(np.shape(M_SC))
    zadd = zz + zzz
    #
    #%subplot(2,2,1), imagesc(zz)
    #%subplot(2,2,2),
    figure()
    title("X")
    imshow(zadd, interpolation='nearest', origin='upper')
    colorbar()
    #%subplot(2,2,3), imagesc(zzz>0)
    #%subplot(2,2,4), imagesc(l+r>0)
    show()
