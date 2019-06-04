#!/usr/bin/env python

import fsnapshot
import sys
from numpy import *
from pylab import *

d = '../prob2'

pf_hi  = d+'/plt98566'
pf_020 = d+'/plt00097'
pf_040 = d+'/plt00386'
pf_080 = d+'/plt01541'
pf_160 = d+'/plt06161'
pf_320 = d+'/plt24642'

def read_data(plotfile):
    (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile)
    t = fsnapshot.fplotfile_get_time(plotfile)

    (xmin, xmax, ymin, ymax, zmin, zmax) = fsnapshot.fplotfile_get_limits(plotfile)

    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    x = xmin + arange( (nx), dtype=float64 )*dx
    y = ymin + arange( (ny), dtype=float64 )*dy

    u = zeros( (nx, ny), dtype=float64)
    (u, err) = fsnapshot.fplotfile_get_data_2d(plotfile, 'phi', u)
    u = transpose(u)

    return t, x, y, u

t_hi, x_hi, y_hi, u_hi = read_data(pf_hi)
t_020, x_020, y_020, u_020 = read_data(pf_020)
t_040, x_040, y_040, u_040 = read_data(pf_040)
t_080, x_080, y_080, u_080 = read_data(pf_080)
t_160, x_160, y_160, u_160 = read_data(pf_160)
t_320, x_320, y_320, u_320 = read_data(pf_320)

nhi = size(x_hi)

utmp = u_hi[0::nhi/20,0::nhi/20] - u_020
err0_020 = abs(utmp).max()
err2_020 = sqrt(sum(utmp**2)*(x_020[1]-x_020[0])*(y_020[1]-y_020[0]))

utmp = u_hi[0::nhi/40,0::nhi/40] - u_040
err0_040 = abs(utmp).max()
err2_040 = sqrt(sum(utmp**2)*(x_040[1]-x_040[0])*(y_040[1]-y_040[0]))

utmp = u_hi[0::nhi/80,0::nhi/80] - u_080
err0_080 = abs(utmp).max()
err2_080 = sqrt(sum(utmp**2)*(x_080[1]-x_080[0])*(y_080[1]-y_080[0]))

utmp = u_hi[0::nhi/160,0::nhi/160] - u_160
err0_160 = abs(utmp).max()
err2_160 = sqrt(sum(utmp**2)*(x_160[1]-x_160[0])*(y_160[1]-y_160[0]))

utmp = u_hi[0::nhi/320,0::nhi/320] - u_320
err0_320 = abs(utmp).max()
err2_320 = sqrt(sum(utmp**2)*(x_320[1]-x_320[0])*(y_320[1]-y_320[0]))

print 'L0 norm errors:', err0_020, err0_040, err0_080, err0_160, err0_320
print 'p0 orders using                ', \
    log(err0_020/err0_040)/log(2.0), \
    log(err0_040/err0_080)/log(2.0), \
    log(err0_080/err0_160)/log(2.0), \
    log(err0_160/err0_320)/log(2.0)

print ''
print 'L2 norm errors:', err2_020, err2_040, err2_080, err2_160, err2_320
print 'p2 orders using                ', \
    log(err2_020/err2_040)/log(2.0), \
    log(err2_040/err2_080)/log(2.0), \
    log(err2_080/err2_160)/log(2.0), \
    log(err2_160/err2_320)/log(2.0)
