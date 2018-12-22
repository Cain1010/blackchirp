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
from scipy.optimize import brentq
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


class BlackChirpMotorScan:
    def __init__(self,number,pretty=False,path=None):
        self.pretty = pretty
        if number < 1:
            raise ValueError('Experiment number must be > 0')
            
        self.d_number = number
        
        millions = number // 1000000
        thousands = number // 1000
        
        if path is None:
            path = "/home/data/blackchirp/experiments/" + str(millions) \
                   + "/" + str(thousands) + "/" + str(number)
        
        fid = open(path+'/'+'%i.mdt'%number,'rb') 
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
        
        self.rawdata = np.zeros((self.Z_size,self.Y_size,self.X_size,self.t_size))
        self.data = np.zeros((self.Z_size,self.Y_size,self.X_size,self.t_size))
        
        
        for i in range(self.Z_size):
            for j in range(self.Y_size):
                for k in range(self.X_size):
                    for l in range(self.t_size):
                        unpacked = struct.unpack(">d", fid.read(struct.calcsize(">d")))[0]
                        self.rawdata[i,j,k,l] = unpacked
                    fid.read(4)
                fid.read(4)
            fid.read(4)
        fid.close()
        self.data = self.rawdata

    def timetrace(self):
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)
        z = 0
        y = 0
        x = 0
        s = self.data[z,y,x,:]
        t = range(len(s))
        l, = plt.plot(t, s, lw=2, color='red')
        ax.set_ylim([0,1])
        
        
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
#            ax.set_ylim([min(d)*0.9,max(d)*1.1])
#            ax.set_ylim([0,1])
        
        sX.on_changed(update)
        sY.on_changed(update)
        sZ.on_changed(update)
    
        
    def pressuretrace(self,label='XY'):
        fig, ax = plt.subplots()
        ax.set_ylim([0,1])
        plt.subplots_adjust(bottom=0.25)
        z = 0
        y = int(self.Y_size/2)-1
        x = int(self.X_size/2)-1
        t = 34
        
        if label == 'XY':
            s = self.data[:,y,x,t]
            t = range(len(s))
            l, = plt.plot(t, s, lw=2, color='red',marker='o')
            ax1 = plt.axes([0.25, 0.14, 0.65, 0.03])#, facecolor=axcolor)
            s1 = Slider(ax1, 'X', 0, self.X_size-1, valinit=x, valfmt='%i')#, valstep=delta_f)
            ax2 = plt.axes([0.25, 0.09, 0.65, 0.03])#, facecolor=axcolor)
            s2 = Slider(ax2, 'Y', 0, self.Y_size-1, valinit=y, valfmt='%i')#, valstep=delta_f)
            def update(val):
                x = int('%i'%s1.val)
                y = int('%i'%s2.val)
                t = int('%i'%st.val)
                d = self.data[:,y,x,t]
                l.set_ydata(d)
                fig.canvas.draw_idle()
#                ax.set_ylim([min(d)*0.9,max(d)*1.1])
        elif label == 'XZ':
            s = self.data[z,:,x,t]
            t = range(len(s))
            l, = plt.plot(t, s, lw=2, color='red',marker='o')
            ax1 = plt.axes([0.25, 0.14, 0.65, 0.03])#, facecolor=axcolor)
            s1 = Slider(ax1, 'X', 0, self.X_size-1, valinit=x, valfmt='%i')#, valstep=delta_f)
            ax2 = plt.axes([0.25, 0.09, 0.65, 0.03])#, facecolor=axcolor)
            s2 = Slider(ax2, 'Z', 0, self.Z_size-1, valinit=z, valfmt='%i')#, valstep=delta_f)
            def update(val):
                x = int('%i'%s1.val)
                z = int('%i'%s2.val)
                t = int('%i'%st.val)
                d = self.data[z,:,x,t]
                l.set_ydata(d)
                fig.canvas.draw_idle()
#                ax.set_ylim([min(d)*0.9,max(d)*1.1])
        elif label == 'YZ':
            s = self.data[z,y,:,t]
            t = range(len(s))
            l, = plt.plot(t, s, lw=2, color='red',marker='o')
            ax1 = plt.axes([0.25, 0.14, 0.65, 0.03])#, facecolor=axcolor)
            s1 = Slider(ax1, 'Y', 0, self.Y_size-1, valinit=y, valfmt='%i')#, valstep=delta_f)
            ax2 = plt.axes([0.25, 0.09, 0.65, 0.03])#, facecolor=axcolor)
            s2 = Slider(ax2, 'Z', 0, self.Z_size-1, valinit=z, valfmt='%i')#, valstep=delta_f)
            def update(val):
                y = int('%i'%s1.val)
                z = int('%i'%s2.val)
                t = int('%i'%st.val)
                d = self.data[z,y,:,t]
                l.set_ydata(d)
                fig.canvas.draw_idle()
#                ax.set_ylim([min(d)*0.9,max(d)*1.1])
        
        
        
        
        axt = plt.axes([0.25, 0.04, 0.65, 0.03])#, facecolor=axcolor)
        st = Slider(axt, 't', 0, self.t_size-1, valinit=34, valfmt='%i')#, valstep=delta_f)
        

        
        s1.on_changed(update)
        s2.on_changed(update)
        st.on_changed(update)
        
#        rax = plt.axes([0.025, 0.05, 0.15, 0.15])#, facecolor=axcolor)
#        radio = RadioButtons(rax, ('XY', 'XZ', 'YZ'), active=0)
#        
#        def axesSwap(label):
#            if label == 'XY':
#                ax1 = plt.axes([0.25, 0.14, 0.65, 0.03])#, facecolor=axcolor)
#                s1 = Slider(ax1, 'X', 0, self.X_size-1, valinit=x, valfmt='%i')
#                ax2 = plt.axes([0.25, 0.09, 0.65, 0.03])#, facecolor=axcolor)
#                s2 = Slider(ax2, 'Y', 0, self.Y_size-1, valinit=y, valfmt='%i')
#            if label == 'XZ':
#                ax1 = plt.axes([0.25, 0.14, 0.65, 0.03])#, facecolor=axcolor)
#                s1 = Slider(ax1, 'X', 0, self.X_size-1, valinit=x, valfmt='%i')
#                ax2 = plt.axes([0.25, 0.09, 0.65, 0.03])#, facecolor=axcolor)
#                s2 = Slider(ax2, 'Z', 0, self.Z_size-1, valinit=z, valfmt='%i')
#            if label == 'YZ':
#                ax1 = plt.axes([0.25, 0.14, 0.65, 0.03])#, facecolor=axcolor)
#                s1 = Slider(ax1, 'Y', 0, self.Y_size-1, valinit=y, valfmt='%i')
#                ax2 = plt.axes([0.25, 0.09, 0.65, 0.03])#, facecolor=axcolor)
#                s2 = Slider(ax2, 'Z', 0, self.Z_size-1, valinit=z, valfmt='%i')
#        radio.on_clicked(axesSwap)    

    def XYCrosssection(self):
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.25)
        z = 0
        t = 34
        CS = self.data[z,:,:,t]
        l = ax.imshow(CS,cmap='hot',vmin=0,vmax=1)
        
        axZ = plt.axes([0.25, 0.1, 0.65, 0.03])#, facecolor=axcolor)
        axt = plt.axes([0.25, 0.15, 0.65, 0.03])#, facecolor=axcolor)
        
        sZ = Slider(axZ, 'Z', 0., self.Z_size-1, valinit=z, valfmt='%i')#, valstep=delta_f)
        st = Slider(axt, 't', 0, self.t_size-1, valinit=t, valfmt='%i')
#        ax.xticks([])
        
        def update(val):
            z = int('%i'%sZ.val)
            t = int('%i'%st.val)
#            print(i,j,k,l)
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
        if self.pretty == True:
            l = ax.imshow(CS,cmap='hot',aspect=.1, interpolation='bilinear',vmin=0,vmax=1)
        else:
            l = ax.imshow(CS,cmap='hot',aspect=.1, interpolation='nearest',vmin=0,vmax=1)
        
        axX = plt.axes([0.17, 0.1, 0.65, 0.03])#, facecolor=axcolor)
        axt = plt.axes([0.17, 0.15, 0.65, 0.03])#, facecolor=axcolor)
        
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
        
    def YZ3DPlot(self,x=9,t=37):
        Y = []
        Z = []
        P = []
        for i in range(self.Z_size):
            for j in range(self.Y_size):
                Y.append(j)
                Z.append(i)
                P.append(self.data[i,j,x,t])
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_wireframe(Y,Z,P)
                           
#    def MachNumbers(self,gamma=1.66):
#        M = np.arange(0.0001,10,.000001)
#        a = (((gamma+1)*M**2)/2)**(gamma/(gamma-1))
#        b = ((gamma+1)/(2*gamma*M**2-(gamma-1)))**(1/(gamma-1)) 
#        Pratio = a*b
#        mtp = {}
#        for i in range(len(M)):
#            mtp[round(Pratio[i],4)] = M[i]
#        self.MachtoPres = mtp
#        self.MachTest = [Pratio,M]
    def getData(self,x,y,z,t):
        pi = self.rawdata[x,y,z,t]
        ps = 0.0333
        pr = round(pi/ps,4)
        return(pr)
        
    def MachNumbers(self,gamma=1.66):
        Pratio = []
        for i in range(self.Z_size):
            for j in range(self.Y_size):
                for k in range(self.X_size):
                    for l in range(self.t_size):
                        Pratio.append(self.getData(i,j,k,l))
        mtp = {}
        for i in range(len(Pratio)):
            PR = Pratio[i]
            try:
                mtp[Pratio[i]] = np.real(brentq((MachNumber),1,200.,args=(Pratio[i],gamma)))
            except TypeError:
                print(PR)
            except ValueError:
                print('outside bounds of equation')
        self.MachtoPres = mtp
        
    def PtoM(self,PR,gamma=1.66):
       try: 
           M = np.real(brentq((MachNumber),0.0000,200.,args=(PR,gamma)))
       except TypeError:
            print(PR)
#        M = brentq((MachNumber),1.,200,args=(0.67,.027,gamma))
       return (M)
    
def MachNumber(M,PR, gamma):
    a = (((gamma+1)*M**2)/2)**(gamma/(gamma-1))
    b = ((gamma+1)/(2*gamma*M**2-(gamma-1)))**(1/(gamma-1))
    return (a*b-(PR))



test = BlackChirpMotorScan(37,True)
#test.YZ3DPlot()
#test.XYCrosssection()
test.timetrace()
#test.YZCrosssection()
#test.pressuretrace()
#test.pressuretrace('XZ')
#test.pressuretrace('YZ')
#test.PtoM()

#test.MachNumbers()

#tmp = test.rawdata[0,9,:,34]
#tmp2 = np.diff(tmp)
#dtmp = []
#for i in range(1,len(tmp)-1):
#    dtmp.append(tmp[i+1]-tmp[i-1])
#    



























