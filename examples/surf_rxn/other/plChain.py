#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.1
#                          Copyright (2021) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
#     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
#     retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is open source software: you can redistribute it and/or modify
#     it under the terms of BSD 3-Clause License
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     BSD 3 Clause License for more details.
#
#     You should have received a copy of the BSD 3 Clause License
#     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
#
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================

import numpy as npy
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
import optparse

global cname, ch, splh, chpl, figh

# command-line options
parser = optparse.OptionParser()
parser.add_option("--cn", dest="cname", default="chain.dat", type="string",
                  help="Name of chain file")
parser.add_option("--nsub", dest="nsub", default="5", type="int",
                  help="No. of subplots")
parser.add_option("--fx", dest="fx", default="6", type="int",
                  help="figure x-size")
parser.add_option("--fy", dest="fy", default="10", type="int",
                  help="figure y-size")
parser.add_option("--skip", dest="nskip", default="5000", type="int",
                  help="skip first # of states")
options, remainder = parser.parse_args()

cname  = options.cname
nsub   = options.nsub
fsizex = options.fx
fsizey = options.fy
nskip = options.nskip

# load chain
ch=npy.genfromtxt(cname)

# no of variables
ist=1;
cend=2;
nvars=ch.shape[1]-ist-cend

# compute number of figures
nfigs=nvars/nsub
if nvars%nsub != 0:
    nfigs += 1;

# create figure handles
figh=[]
for k in range(nfigs):
    figh.append(figure(figsize=(fsizex,fsizey)))

# create subplot handles
splh=[]
for i in range(nvars):
    ifg = i/nsub;
    jsb = (i%nsub)+1
    splh.append(figh[ifg].add_subplot(nsub,1,jsb))

# plot chain
chpl=[]
for i in range(nvars):
    pl1,=splh[i].plot(ch[nskip:,0],ch[nskip:,i+1])
    chpl.append(pl1)
    splh[i].set_ylabel("v"+str(i+1))

# define update action
def reset(event):
    global cname,ch,splh,chpl,figh
    ch=npy.genfromtxt(cname)
    for i in range(nvars):
        chpl[i].set_xdata(ch[nskip:,0])
        chpl[i].set_ydata(ch[nskip:,i+1])
        splh[i].set_xlim([ch[nskip:,0].min(),ch[nskip:,0].max()])
        splh[i].set_ylim([ch[nskip:,i+1].min(),ch[nskip:,i+1].max()])
    for fig in figh:
        fig.canvas.draw()
    print "plots were updated"

figb=figure(figsize=(1,0.5))
resetax = axes([0.1, 0.2, 0.8, 0.6])
axcolor = 'lightgoldenrodyellow'
button = Button(resetax, 'Load chain', color=axcolor, hovercolor='0.975')
button.on_clicked(reset)
show()

# if plw == "chn":
#     for i in range(1,nvars+1):
#         fig=plt.figure()
#         plt.plot(ch[chst::,i])
# elif plw == "hist":
#     for i in range(1,nvars+1):
#         fig=plt.figure()
#         a11=plt.hist(ch[chst::,i],nbins)
#         plt.savefig("hist_"+str(i+1)+".eps")
# else:
#     for i in range(1,nvars+1):
#         fig=plt.figure()
#         plt.plot(ch[chst::,i],ch[chst::,nvars+2],linestyle='none',marker='.')


# # define update action
# def reset(event):
#     global cname,ch,splh,chpl
#     ch=npy.genfromtxt(cname)
#     print 'button was pushed ',ch[:,0].max()

# figb=figure(figsize=(1,0.5))
# resetax = axes([0.1, 0.2, 0.8, 0.6])
# axcolor = 'lightgoldenrodyellow'
# button = Button(resetax, 'Load chain', color=axcolor, hovercolor='0.975')
# button.on_clicked(reset)
# show()

