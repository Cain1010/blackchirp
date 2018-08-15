# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 17:47:49 2018

@author: zsbuchanan
"""

import struct
import numpy as np
import scipy.signal as spsig
import shelve
import os
import bcfitting as bcfit
import concurrent.futures
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons


class BlackChirpMotorScan:
    def __init__(self,number,pretty=False,path=None):
        
        if number < 1:
            raise ValueError('Experiment number must be > 0')
            
        self.d_number = number
        
        millions = number // 1000000
        thousands = number // 1000
        
        if path is None:
            path = "/home/data/blackchirp/experiments/" + str(millions) \
                   + "/" + str(thousands) + "/" + str(number)
        
        fid = open(path+'/%i.mdt'%number,'rb') 
        buffer = fid.read(4)
        ms_len = struct.unpack(">I",buffer)
        buffer = fid.read(ms_len[0])
        magic_string = buffer.decode('ascii')
        buffer = fid.read(4)
        self.Z_size = struct.unpack(">I", buffer)[0]
        buffer = fid.read(4)
        self.Y_size = struct.unpack(">I", buffer)[0]
        buffer = fid.read(4)
        self.X_size = struct.unpack(">I", buffer)[0]
        buffer = fid.read(4)
        self.t_size = struct.unpack(">I", buffer)[0]
        
        self.data = np.zeros((self.Z_size,self.Y_size,self.X_size,self.t_size))
        
        
        for i in range(self.Z_size):
            for j in range(self.Y_size):
                for k in range(self.X_size):
                    for l in range(self.t_size):
                        unpacked = struct.unpack(">d", fid.read(struct.calcsize(">d")))[0]
                        self.data[i,j,k,l] = unpacked
                    fid.read(4)
                fid.read(4)
            fid.read(4)
        fid.close()
        
        fid = open(path+'%i.hdr'%number)

    def timetrace(self):
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)
        z = 0
        y = 0
        x = 0
        s = self.data[z,y,x,:]
        t = range(len(s))
        l, = plt.plot(t, s, lw=2, color='red')
        
        
#        axcolor = 'lightgoldenrodyellow'
        axX = plt.axes([0.25, 0.14, 0.65, 0.03])#, facecolor=axcolor)
        sX = Slider(axX, 'X', 0, self.X_size-1, valinit=x, valfmt='%i')#, valstep=delta_f)
        axY = plt.axes([0.25, 0.09, 0.65, 0.03])#, facecolor=axcolor)
        sY = Slider(axY, 'Y', 0, self.Y_size-1, valinit=x, valfmt='%i')#, valstep=delta_f)
        axZ = plt.axes([0.25, 0.04, 0.65, 0.03])#, facecolor=axcolor)
        sZ = Slider(axZ, 'Z', 0, self.Z_size-1, valinit=x, valfmt='%i')#, valstep=delta_f)
        
        def update(val):
            x = int('%i'%sX.val)
            y = int('%i'%sY.val)
            z = int('%i'%sZ.val)
            d = self.data[z,y,x,:]
            l.set_ydata(d)
            fig.canvas.draw_idle()
            ax.set_ylim([min(d)*0.9,max(d)*1.1])
        
        sX.on_changed(update)
        sY.on_changed(update)
        sZ.on_changed(update)

    def XYCrosssection(self):
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)
        z = 0
        t = 34
        CS = self.data[z,:,:,t]
        l = ax.imshow(CS,cmap='hot')
        
        axZ = plt.axes([0.25, 0.1, 0.65, 0.03])#, facecolor=axcolor)
        axt = plt.axes([0.25, 0.15, 0.65, 0.03])#, facecolor=axcolor)
        
        sZ = Slider(axZ, 'Z', 0., self.Z_size-1, valinit=z, valfmt='%i')#, valstep=delta_f)
        st = Slider(axt, 't', 0, self.t_size-1, valinit=t, valfmt='%i')
        
        
        def update(val):
            z = int('%i'%sZ.val)
            t = int('%i'%st.val)
            l.set_data(self.data[z,:,:,t])
#            print('updated to: ', z, t)
            fig.canvas.draw_idle()
            
        sZ.on_changed(update)
        st.on_changed(update)
        
    def YZCrosssection(self):
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)
        x = 0
        t = 34
        CS = np.transpose(self.data[:,:,x,t])
        l = ax.imshow(CS,cmap='hot',aspect=.2, interpolation='nearest')
        
        axX = plt.axes([0.25, 0.1, 0.65, 0.03])#, facecolor=axcolor)
        axt = plt.axes([0.25, 0.15, 0.65, 0.03])#, facecolor=axcolor)
        
        sX = Slider(axX, 'X', 0., self.X_size-1, valinit=x, valfmt='%i')#, valstep=delta_f)
        st = Slider(axt, 't', 0, self.t_size-1, valinit=t, valfmt='%i')
        
        
        def update(val):
            x = int('%i'%sX.val)
            t = int('%i'%st.val)
            l.set_data(np.transpose(self.data[:,:,x,t]))
#            print('updated to: ', z, t)
            fig.canvas.draw_idle()
            
        sX.on_changed(update)
        st.on_changed(update)


test = BlackChirpMotorScan(36)
#test.XYCrosssection()
#test.timetrace()
test.YZCrosssection()
