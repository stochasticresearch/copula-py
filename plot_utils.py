#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#******************************************************************************
#* 
#* Copyright (C) 2015  Kiran Karra <kiran.karra@gmail.com>
#*
#* This program is free software: you can redistribute it and/or modify
#* it under the terms of the GNU General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or
#* (at your option) any later version.
#*
#* This program is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#* GNU General Public License for more details.
#*
#* You should have received a copy of the GNU General Public License
#* along with this program.  If not, see <http://www.gnu.org/licenses/>.
#******************************************************************************

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

def plot_3d(X,Y,Z, titleStr):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.xlabel('U1')
    plt.ylabel('U2')
    plt.title(titleStr)
    plt.show()
    
def pairs(X, titleStr):
    # TODO: make this generic, rightnow specific to 3 variate data
    f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3)    
    ax2.scatter(X[:,0],X[:,1])
    ax2.set_title(titleStr + ' 1 : 2')
    ax2.grid()
    
    ax3.scatter(X[:,0],X[:,2])
    ax3.set_title(titleStr + ' 1 : 3')
    ax3.grid()
    
    ax4.scatter(X[:,1],X[:,0])
    ax4.set_title(titleStr + ' 2 : 1')
    ax4.grid()
    
    ax6.scatter(X[:,1],X[:,2])
    ax6.set_title(titleStr + ' 2 : 3')
    ax6.grid()
    
    ax7.scatter(X[:,2],X[:,0])
    ax7.set_title(titleStr + ' 3 : 1')
    ax7.grid()
    
    ax8.scatter(X[:,2],X[:,1])
    ax8.set_title(titleStr + ' 3 : 2')
    ax8.grid()
    
    plt.show()