'''
 Make a vortex ready for laser cutting. 
 

'''

from mpl_toolkits import mplot3d
import matplotlib
import numpy as np
from matplotlib import pyplot as plt


# Implementation of Z-buffer Algorithm
# 4.14.2002.
# Jiwon Hahn

# def pixel(): return 100  # 32x32 grid of pixels
# def inf(): return 'inf' 

# def z_buffer_algo():  
#     zbuff=[[inf()]*pixel() for i in range(pixel())] #initialize z buffer to inf
#     intensity=[[0]*pixel() for i in range(pixel())] #initialize intensity to zero

#     print ("Polygon 1(rectangular):")
#     for x in range(pixel()):    
#         for y in range(pixel()): 
#             if (x>=10 and x<=25 and y>=5 and y<=25): #point inside the rectangular
#                 z_depth=10
#                 if z_depth < zbuff[x][y] or zbuff[x][y]=='inf':
#                     intensity[x][y]=1
#                     zbuff[x][y]=z_depth
#     print ("------------------------Z buffer-------------------------")
#     for x in range(pixel()-1,-1,-1):
#         print (zbuff[x])

#     print ("------------------------Intensity------------------------")
#     for x in range(pixel()-1,-1,-1):
#         print (intensity[x])
#     print ("---------------------------------------------------------")

#     print ("Polygon 2(triangle):")
#     for x in range(pixel()):
#         for y in range(pixel()):
#             if (y<=x and y<=100-3*x and y>= 20-1.0/3*x):
#                 z_depth= 30-3.0/4*x-1.0/4*y
#                 if z_depth < zbuff[x][y] or zbuff[x][y]=='inf':
#                     intensity[x][y]=2
#                     zbuff[x][y]=z_depth
#     print ("------------------------Z buffer-------------------------")
#     for x in range(pixel()-1,-1,-1):
#         print (zbuff[x])

#     print ("------------------------Intensity------------------------")
#     for x in range(pixel()-1,-1,-1):
#         print (intensity[x])
#     print ("---------------------------------------------------------")
                

# # test
# z_buffer_algo()
            

n = 100
theta 		= np.linspace(np.pi/4, 1.0*np.pi, n)
phi 		= np.linspace(0, 2*np.pi, n)
theta, phi 	= np.meshgrid(theta, phi)
c, a = 2.5, 2.0
b, d = 3, 1.9 
y = (c + a*np.cos(theta)) * np.cos(phi)
z = (b + d*np.cos(theta)) * np.sin(phi)
x = a * np.sin(theta)

fig = plt.figure()
ax2 = fig.add_subplot(111,projection='3d')
# ax1 = fig.add_subplot(121, projection='3d')
# ax1.set_zlim(-3,3)
# ax1.contour3D(x,y,z, 50,cmap='binary')
# # ax1.plot_surface(x, y, z, rstride=5, cstride=5, color='k', edgecolors='w')
# ax1.view_init(30, 137)
# ax2 = fig.add_subplot(122, projection='3d')
ax2.set_zlim(-3,3)

ax2.plot_surface(x, y, z, rstride=3, cstride=8, color='white', edgecolors='k',shade=False)

# surf = ax2.plot_surface(x,y,z, rstride=2, cstride=2, shade=False, cmap="jet", linewidth=1)
# surf.set_edgecolors(surf.to_rgba(surf._A))
# surf.set_facecolors("white")

ax2.view_init(350,20)
ax2.set_axis_off()
plt.savefig("test.svg", format="svg")
plt.show()

# ele = 55
# azm = 55
# dst = 10 
# ax2.view_init(elev=ele, azim=azm) #Works!
# ax2.dist=dst 
# 
# elev=10., azim=ii
# rotate the axes and update
# for angle in range(0, 360,):
#     # ax2.view_init(55, angle)
#     ax2.view_init(350, angle)     # 180+40, elev 330-360 azim 20
#     plt.draw()
#     print(angle)
#     plt.pause(.1)
# 


# import bpy
# from bpy import context
# filepath = "test.obj"
# python wavefront object? 
# with open(filepath, 'w') as f:
#     f.write("# OBJ file\n")
#     for v in mesh:
#         f.write("v %.4f %.4f %.4f\n" % (x[v],y[v],z[v]))
#     for p in mesh.polygons:
#         f.write("f")
#         for i in p.vertices:
#             f.write(" %d" % (i + 1))
#         f.write("\n")
# import os



