#!/usr/bin/env python

import fsnapshot
import sys
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import mpl_toolkits.axes_grid as mpl_ag

pf_wide = '../prob3_chain/plt00386'
pf_narr = '../prob3/plt00386'

#pf_wide = '../prob3_chain/plt01541'
#pf_narr = '../prob3/plt01541'

# pf_wide

(nx, ny, nz) = fsnapshot.fplotfile_get_size(pf_wide)
t = fsnapshot.fplotfile_get_time(pf_wide)

(xmin, xmax, ymin, ymax, zmin, zmax) = fsnapshot.fplotfile_get_limits(pf_wide)

dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
x = xmin + arange( (nx), dtype=float64 )*dx
y = ymin + arange( (ny), dtype=float64 )*dy

u = zeros( (nx, ny), dtype=float64)
(u, err) = fsnapshot.fplotfile_get_data_2d(pf_wide, 'phi', u)
u = transpose(u)
u_wide = zeros( (nx+1, ny+1), dtype=float64)
u_wide[0:nx,0:ny] = u
u_wide[nx,0:ny] = u[0,0:ny]
u_wide[0:nx,ny] = u[0:nx,0]
u_wide[nx,ny] = u[0,0]

# pf_narr

(nx, ny, nz) = fsnapshot.fplotfile_get_size(pf_narr)
t = fsnapshot.fplotfile_get_time(pf_narr)

(xmin, xmax, ymin, ymax, zmin, zmax) = fsnapshot.fplotfile_get_limits(pf_narr)

dx = (xmax - xmin) / nx
dy = (ymax - ymin) / ny
x = xmin + arange( (nx), dtype=float64 )*dx
y = ymin + arange( (ny), dtype=float64 )*dy

(u, err) = fsnapshot.fplotfile_get_data_2d(pf_narr, 'phi', u)
u = transpose(u)
u_narr = zeros( (nx+1, ny+1), dtype=float64)
u_narr[0:nx,0:ny] = u
u_narr[nx,0:ny] = u[0,0:ny]
u_narr[0:nx,ny] = u[0:nx,0]
u_narr[nx,ny] = u[0,0]

#####################################

fig = figure(2, figsize=(8,6))

ax1 = fig.add_subplot(121, projection='3d')
mx, my = mgrid[0:xmax+dx:dx, 0:ymax+dy:dy]
surf1 = ax1.plot_surface(mx, my, u_wide, cmap='jet', rstride=1, cstride=1)
ax1.view_init(15,250)
ax1.set_xlim(0,2.*pi)
ax1.set_ylim(0,2.*pi)
ax1.set_title('Wide Stencil')

ax2 = fig.add_subplot(122, projection='3d')
surf2 = ax2.plot_surface(mx, my, u_narr, cmap='jet', rstride=1, cstride=1)
ax2.view_init(15,250)
ax2.set_xlim(0,2.*pi)
ax2.set_ylim(0,2.*pi)
ax2.set_title('Narrow Stencil')

draw()
show()
