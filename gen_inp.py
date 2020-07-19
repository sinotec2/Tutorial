#!/opt/anaconda3/envs/py27/bin/python
import numpy as np
from pandas import *
import twd97, utm
from scipy.interpolate import griddata
import matplotlib._cntr as cntr
from libtiff import TIFF,libtiff_ctypes
import bisect
import sys,os

#  conda install matplotlib==2.1.0


def getarg():
  """ read the setting of plot from argument(std input)
  -f --fname x,y,c data file name
  -d --dictr(optional) distric level (D5/D6/ALL)
  the resultant file is FNAME.kml, can be plot in google map interface
  """
  import argparse
  ap = argparse.ArgumentParser()
  ap.add_argument("-f", "--fname", required=True, type=str, help="isc_plotfile")
  ap.add_argument("-d", "--dictr", required=False, type=str, help="region of plot")
  args = vars(ap.parse_args())
  return args['fname'], args['dictr']


args = getarg()
fname = args[0]
distr = 'ALL'
xmin,nx,dx,ymin,ny,dy=(float(i) for i in args[1].split())
nx,ny=int(nx),int(ny)
xmax=xmin+(nx-1)*dx
ymax=ymin+(ny-1)*dy

# read the taiwan.tif file, must be in TWD97-m system, 
libtiff_ctypes.suppress_warnings()
cap = TIFF.open('/Users/cybee/taiwan.tif')
d = cap.read_image()
mn=d.shape

# the following parameters are read from linux command: "tiffinfo taiwan.tif"
x0,y0=42033.952431,2878094.424011
dxm,dym=28.884201,-28.884201

# 1-d X/Y coordinates
x1d = np.array([x0+dxm*i for i in range(mn[1])])
y1d = np.array([y0+dym*i for i in range(mn[0])])
y1dr=y1d
y1dr.sort() #resort are dangerous, may affect previous array (y1d)
y1d = np.array([y0+dym*i for i in range(mn[0])])

#extract by given domain defined by -d "..."
I1=bisect.bisect_left(x1d,xmin)-2
I2=bisect.bisect_left(x1d,xmax)+2
J1=bisect.bisect_left(y1dr,ymin)-2
J2=bisect.bisect_left(y1dr,ymax)+2
c=d[-J2:-J1,I1:I2].flatten()
if np.min(c)==np.max(c):
  sys.exit('flat terrain chosen!')
x=np.array([x1d[I1:I2] for j in range(J2-J1)]).flatten()
y=np.array([[j for i in range(I2-I1)] for j in y1d[-J2:-J1]]).flatten()

#the domain of meshes must smaller than data domain to avoid extra_polation
x_mesh = np.linspace(xmin, xmax, nx)
y_mesh = np.linspace(ymin, ymax, ny)


# 2-d mesh coordinates, both in TWD97 and WGS84
x_g, y_g = np.meshgrid(x_mesh, y_mesh)

#lat,lon pairs are used in KML locations
ll = np.array([[twd97.towgs84(i, j) for i, j in zip(x_g[k], y_g[k])] for k in xrange(ny)])
lat, lon = (ll[:, :, i] for i in [0, 1])
e, w, s, n = np.max(lon), np.min(lon), np.min(lat), np.max(lat)

# interpolation vector c into meshes
points = np.array(zip(x, y)).astype(float)
# cubic spline may smooth too much,be careful.
#grid_z2 = griddata(points, c, (x_g, y_g))#, method='cubic')
grid_z2 = griddata(points, c, (x_g, y_g), method='linear')
grid_z2 = np.clip(grid_z2,0.,np.max(grid_z2))
# levels size,>10 too thick, <5 too thin
N = 10
levels = np.linspace(0, np.max(grid_z2), N)
col = '#00FF0A #3FFF0A #7FFF0A #BFFF0A #FFFF0A #FECC0A #FD990A #FC660A #FB330A #FA000A'.replace('#', '').split()
if len(col) != N: print 'color scale not right, please redo from http://www.zonums.com/online/color_ramp/'
aa = '28'  # ''28'~ 40%, '4d' about 75%
rr, gg, bb = ([i[j:j + 2] for i in col] for j in [0, 2, 4])
col = [aa + b + g + r for b, g, r in zip(bb, gg, rr)]

# round the values of levels to 1 significant number at least, -2 at least 2 digits
M = int(np.log10(levels[1])) - 1
levels = [round(lev, -M) for lev in levels]

#the Cntr method is valid only in previous version of matplotlib
c = cntr.Cntr(lon, lat, grid_z2)
# the tolerance to determine points are connected to the boundaries
tol = 1E-3
col0 = '4d6ecdcf'
col_line0 = 'cc2d3939'


#writing the KML, see the KML official website
head1 = '<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://earth.google.com/kml/2.2"><Document><name><![CDATA[' + fname + ']]></name>'
st_head = ''
st_med = '</color><width>1</width></LineStyle><PolyStyle><color>'
st_tail = '</color></PolyStyle></Style>'
for i in xrange(N):
  st_head += '<Style id="level' + str(i) + '"><LineStyle><color>' + col[i] + st_med + col[i] + st_tail
head2 = '</styleUrl><Polygon><outerBoundaryIs><LinearRing><tessellate>1</tessellate><coordinates>'
tail2 = '</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>'
line = [head1 + st_head]
# repeat for the level lines
for level in levels[:]:
  nlist = c.trace(level, level, 0)
  segs = nlist[:len(nlist) // 2]
  i = levels.index(level)
  for seg in segs:
    line.append('<Placemark><name>level:' + str(level) + '</name><styleUrl>#level' + str(i) + head2)
    leng = -9999
    for j in xrange(len(seg[:, 0])):
      line.append(str(seg[j, 0]) + ',' + str(seg[j, 1]) + ',0 ')
      if j > 0:
        leng = max(leng, np.sqrt((seg[j, 0] - seg[j - 1, 0]) ** 2 + (seg[j, 1] - seg[j - 1, 1]) ** 2))
    leng0 = np.sqrt((seg[j, 0] - seg[0, 0]) ** 2 + (seg[j, 1] - seg[0, 1]) ** 2)
    ewsn = np.zeros(shape=(4, 2))
    j = -1
    # end points not closed, add coner point(s) to close the polygons.
    if leng0 > leng and leng0 / leng > 5:
      if abs(seg[j, 0] - e) < tol: ewsn[0, 1] = 1
      if abs(seg[0, 0] - e) < tol: ewsn[0, 0] = 1
      if abs(seg[j, 0] - w) < tol: ewsn[1, 1] = 1
      if abs(seg[0, 0] - w) < tol: ewsn[1, 0] = 1
      if abs(seg[j, 1] - s) < tol: ewsn[2, 1] = 1
      if abs(seg[0, 1] - s) < tol: ewsn[2, 0] = 1
      if abs(seg[j, 1] - n) < tol: ewsn[3, 1] = 1
      if abs(seg[0, 1] - n) < tol: ewsn[3, 0] = 1
      if sum(ewsn[1, :] + ewsn[2, :]) == 2: line.append(str(np.min(lon)) + ',' + str(np.min(lat)) + ',0 ')
      if sum(ewsn[1, :] + ewsn[3, :]) == 2: line.append(str(np.min(lon)) + ',' + str(np.max(lat)) + ',0 ')
      if sum(ewsn[0, :] + ewsn[3, :]) == 2: line.append(str(np.max(lon)) + ',' + str(np.max(lat)) + ',0 ')
      if sum(ewsn[0, :] + ewsn[2, :]) == 2: line.append(str(np.max(lon)) + ',' + str(np.min(lat)) + ',0 ')
    # TODO: when contour pass half of the domain,must add two edge points.
    line.append(tail2)
line.append('</Document></kml>')
with open(fname + '.kml', 'w') as f:
  [f.write(i) for i in line]

xy = np.array([[(i, j) for i, j in zip(x_g[k], y_g[k])] for k in xrange(ny)])
with open(fname + '_re.dat','w') as f:
  f.write('RE ELEVUNIT METERS\n')
  for j in range(ny):
    for i in range(nx):
      f.write('RE DISCCART '+str(xy[j,i,0])+' '+str(xy[j,i,1])+' '+str(grid_z2[j,i])+'\n')

#terrain grid file
with open(fname + '_TG.txt','w') as f:
  f.write(str(nx)+' '+str(ny)+' '+str(xmin)+' '+str(xmax)+' '+str(ymin)+' '+str(ymax)+' '+str(dx)+' '+str(dy)+'\n')
  for j in range(ny):
    ele=[str(int(grid_z2[j,i])) for i in range(nx)]
    s=ele[0]
    for i in range(1,nx):
      s+=' '+ele[i] 
    f.write(s+'\n')

#download from eio
llmin=twd97.towgs84(xmin-2000.*dx/100, ymin-2000*dx/100.)
llmax=twd97.towgs84(xmax+2000.*dy/100, ymax+2000*dy/100.)
smax=str(llmax[1])+' '+str(llmax[0])
smin=str(llmin[1])+' '+str(llmin[0])+' '+smax
llSE=str(llmax[1])+' '+str(llmin[0])
llNE=str(llmin[1])+' '+str(llmax[0])+' '+llSE
TIF,DEM,NUL=fname+'.tiff',fname+'.dem','>&/dev/null'
os.system('/opt/anaconda3/bin/eio clip -o '+TIF+' --bounds '+smin+NUL)
os.system('/opt/anaconda3/envs/ncl_stable/bin/gdal_translate -of USGSDEM -ot Float32 -projwin '+llNE+' '+TIF+' '+DEM+NUL)

#generate aermap.inp, 
#domain of rectangular region, dx*3 is due to the resolution of DEM file may not sufficient
llmin=twd97.towgs84(xmin-dx*5, ymin-dx*5)
llmax=twd97.towgs84(xmax+dx*5, ymax+dy*5)
uxmn1,uymn1=utm.from_latlon(llmin[0],llmin[1])[0:2]
uxmn2,uymx1=utm.from_latlon(llmax[0],llmin[1])[0:2]
uxmx1,uymn2=utm.from_latlon(llmin[0],llmax[1])[0:2]
uxmx2,uymx2=utm.from_latlon(llmax[0],llmax[1])[0:2]
co,an,z='   DOMAINXY  ','   ANCHORXY  ',' 51 '
UTMrange=co+str(int(max(uxmn1,uxmn2)))+' '+str(int(max(uymn1,uymn2)))+z+str(int(min(uxmx1,uxmx2)))+' '+str(int(min(uymx1,uymx2)))+z
xmid,ymid=(xmin+xmax)/2., (ymin+ymax)/2.
llanc=twd97.towgs84(xmid,ymid)
uxanc,uyanc=utm.from_latlon(llanc[0],llanc[1])[0:2]
UTMancha=an+str(int(xmid))+' '+str(int(ymid))+' '+str(int(uxanc))+' '+str(int(uyanc))+z+'0'

#change the contain of aermap
text_file = open("aermap.inp", "r")
d=[line for line in text_file]
keywd=[i[3:11] for i in d]
ifile=keywd.index('DATAFILE')
idmxy=keywd.index('DOMAINXY')
ianxy=keywd.index('ANCHORXY')
ihead=keywd.index('ELEVUNIT')
iend=d.index('RE FINISHED\n')

text_file = open("aermap.inp", "w")

x0,y0=xmin,ymin
s=[]
for j in xrange(ny):
  dyj=dy*(float(j)+0.5)
  for i in xrange(nx):
    dxi=dx*(float(i)+0.5)
    s.append('   DISCCART  '+str(round(x0+dxi,1))+' '+str(round(y0+dyj,1))+'\n')
for l in xrange(ihead+1):
  if l == ifile:
    text_file.write( "%s" % '   DATAFILE  '+DEM+'\n')
  elif l == idmxy:
    text_file.write( "%s" % UTMrange+'\n')
  elif l == ianxy:
    text_file.write( "%s" % UTMancha+'\n')
  else:	
    text_file.write( "%s" % d[l])
for l in xrange(len(s)):
    text_file.write( "%s" % s[l])
for l in xrange(iend,len(d)):
    text_file.write( "%s" % d[l])
text_file.close()

