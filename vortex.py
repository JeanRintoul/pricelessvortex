#!/usr/bin/env python
'''Make a vortex ready for laser cutting.'''

# Fix backend (https://stackoverflow.com/questions/31373163)
import matplotlib
matplotlib.use('TkAgg')

from mpl_toolkits import mplot3d
import matplotlib
import numpy as np
import scipy.linalg
from matplotlib import pyplot as plt
from matplotlib import collections


# Implementation of Z-buffer Algorithm
# 4.14.2002.
# Jiwon Hahn

# def pixel(): return 100  # 32x32 grid of pixels
# def inf(): return 'inf'

def z_buffer_algo():
    zbuff=[[inf()]*pixel() for i in range(pixel())] #initialize z buffer to inf
    intensity=[[0]*pixel() for i in range(pixel())] #initialize intensity to zero

    print("Polygon 1(rectangular):")
    for x in range(pixel()):
        for y in range(pixel()):
            if (x>=10 and x<=25 and y>=5 and y<=25): #point inside the rectangular
                z_depth=10
                if z_depth < zbuff[x][y] or zbuff[x][y]=='inf':
                    intensity[x][y]=1
                    zbuff[x][y]=z_depth
    print("------------------------Z buffer-------------------------")
    for x in range(pixel()-1,-1,-1):
        print(zbuff[x])

    print("------------------------Intensity------------------------")
    for x in range(pixel()-1,-1,-1):
        print(intensity[x])
    print("---------------------------------------------------------")

    print("Polygon 2(triangle):")
    for x in range(pixel()):
        for y in range(pixel()):
            if (y<=x and y<=100-3*x and y>= 20-1.0/3*x):
                z_depth= 30-3.0/4*x-1.0/4*y
                if z_depth < zbuff[x][y] or zbuff[x][y]=='inf':
                    intensity[x][y]=2
                    zbuff[x][y]=z_depth
    print("------------------------Z buffer-------------------------")
    for x in range(pixel()-1,-1,-1):
        print(zbuff[x])

    print("------------------------Intensity------------------------")
    for x in range(pixel()-1,-1,-1):
        print(intensity[x])
    print("---------------------------------------------------------")


# # test
# z_buffer_algo()

def surface(n):
    theta = np.linspace(np.pi/4, 1.0*np.pi, n)
    phi = np.linspace(0, 2*np.pi, n)
    theta, phi = np.meshgrid(theta, phi)
    c, a = 2.5, 2.0
    b, d = 3, 1.9
    y = (c + a*np.cos(theta)) * np.cos(phi)
    z = (b + d*np.cos(theta)) * np.sin(phi)
    x = a * np.sin(theta)
    return np.stack([x, y, z], axis=-1)


def mesh(X, strides=(1,1)):
    assert X.shape[-1] == 3
    s0, s1 = strides
    n0, n1, _ = X.shape
    X = X.reshape(n0*n1, 3)
    i0 = np.arange(n0-1)[:,None]
    i1 = np.arange(n1-1)[None,:]
    i = i0*n1 + i1
    tris = np.concatenate([np.stack([i+n1,i,i+1],axis=-1),
                           np.stack([i+1,i+1+n1,i+n1],axis=-1)]).reshape(-1,3)
    segs = np.concatenate([np.stack([i[::s0],i[::s0]+1], axis=-1).reshape(-1,2),
                           np.stack([i[:,::s1],i[:,::s1]+n1], axis=-1).reshape(-1,2)])
    return X, tris, segs


def clip_hidden(tris, segs):
    # Require double precision
    assert tris.dtype == np.float64
    assert segs.dtype == np.float64

    # Randomly rotate by epsilon to avoid axis degeneracy
    A = 1e-4*np.random.randn(3,3)
    A = A - A.T
    Q = scipy.linalg.expm(A)
    tris = np.dot(tris, Q.T)
    segs = np.dot(segs, Q.T)

    # Ensure triangles are positively oriented
    def normals(tris):
      X0, X1, X2 = tris.swapaxes(0, 1)
      return np.cross(X1-X0, X2-X0)
    tris = tris[normals(tris)[:,2] != 0]
    inverted, = np.nonzero(normals(tris)[:,2] < 0)
    tris[inverted] = tris[inverted][:, ::-1]
    normals = normals(tris)
    assert np.all(normals[:,2] > 0)
    normals /= np.linalg.norm(normals, ord=2, axis=-1, keepdims=True)

    # Clip each segment
    clipped = []
    edges = np.asarray([[0,1],[1,2],[2,0]])
    for S in segs:
        center = (S[0] + S[1]) / 2
        S = S[:,:2]
        DS = S[1] - S[0]

        # At this point, interpolation range [0,1] is visible on the segment.
        # We model this as two events going from count 0 to 1 at time 0, and count 1 to 0 at time 1.
        times = [[0, 1]]
        steps = [[1, -1]]

        # Add events for all triangles
        front = ((center - tris[:,0]) * normals).sum(axis=-1) < -1e-6
        T = tris[front,:,:2]
        E = T[np.arange(len(T))[:,None,None], edges]
        DX = E[:,:,1] - E[:,:,0]
        det = np.cross(DX, DS)
        s = np.cross(DX, E[:,:,0] - S[0]) / (det + 1e-30)
        tlo = np.where(det > 0, s, -np.inf).max(axis=-1)
        thi = np.where(det < 0, s, np.inf).min(axis=-1)
        i, = np.nonzero(tlo < thi)
        times.append(np.stack([tlo[i], thi[i]], axis=-1).ravel())
        steps.append(np.tile([-1, 1], len(i)))

        # Sort events by time
        times = np.concatenate(times)
        steps = np.concatenate(steps)
        p = np.argsort(times)
        times = times[p]
        steps = steps[p]
        values = np.cumsum(steps)

        # The clipped segments are all intervals where the total step is 1
        starts, = np.nonzero(np.logical_and(values == 1, steps == 1))
        stops, = np.nonzero(np.logical_and(values == 0, steps == -1))
        s = np.stack([times[starts], times[stops]], axis=-1)
        clipped.extend(S[0] + DS*s[...,None])
    clipped = np.asarray(clipped)
    return clipped


def plot_surface(X, return_proj=False):
    fig = plt.figure()
    ax2 = fig.add_subplot(111,projection='3d')
    # ax1 = fig.add_subplot(121, projection='3d')
    # ax1.set_zlim(-3,3)
    # ax1.contour3D(x,y,z, 50,cmap='binary')
    # # ax1.plot_surface(x, y, z, rstride=5, cstride=5, color='k', edgecolors='w')
    # ax1.view_init(30, 137)
    # ax2 = fig.add_subplot(122, projection='3d')
    ax2.set_zlim(-3,3)

    ax2.plot_surface(*X.T, rstride=3, cstride=8, color='white', edgecolors='k',shade=False)

    # surf = ax2.plot_surface(x,y,z, rstride=2, cstride=2, shade=False, cmap="jet", linewidth=1)
    # surf.set_edgecolors(surf.to_rgba(surf._A))
    # surf.set_facecolors("white")

    ax2.view_init(350,20)
    ax2.set_axis_off()
    if return_proj:
        return ax2.get_proj()
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


def plot_segs(Xsegs):
    lines = collections.LineCollection(Xsegs[...,:2])
    fig, axis = plt.subplots(1)
    axis.add_collection(lines)
    axis.autoscale()
    axis.margins(0.1)
    axis.set_axis_off()
    print('Show!')
    plt.show()


def write_obj(X):
    import bpy
    from bpy import context
    filepath = "test.obj"
    # python wavefront object?
    with open(filepath, 'w') as f:
        f.write("# OBJ file\n")
        for v in mesh:
            f.write("v %.4f %.4f %.4f\n" % (x[v],y[v],z[v]))
        for p in mesh.polygons:
            f.write("f")
            for i in p.vertices:
                f.write(" %d" % (i + 1))
            f.write("\n")
    import os


def transform(P, X):
    X = np.concatenate([X, np.ones((len(X),1))], axis=-1)
    X = X @ P.T
    X = X[:,:3] / X[:,3:]
    X[:,2] *= -1
    return X


if __name__ == '__main__':
    cmd = 'clip'
    if cmd == 'surface':
        X = surface(n=96)
        plot_surface(X)
    elif cmd == 'raw':
        X = surface(n=96)
        P = plot_surface(X, return_proj=True)
        X, tris, segs = mesh(X, strides=(8,3))
        X = transform(P, X)
        plot_segs(X[segs])
    elif cmd == 'clip':
        X = surface(n=96)
        P = plot_surface(X, return_proj=True)
        X, tris, segs = mesh(X, strides=(8,3))
        X = transform(P, X)
        clipped = clip_hidden(X[tris], X[segs])
        plot_segs(clipped)
    elif cmd == 'single':
        X = np.asarray([[0,1,0],[1,-1,0],[-1,-1,0]], dtype=np.float64)
        tris = np.asarray([[0,1,2]])
        segs = np.asarray([[[-2,y,y],[2,y,y]] for y in np.linspace(-1,1)], dtype=np.float64)
        clipped = clip_hidden(X[tris], segs)
        plot_segs(clipped)
    elif cmd == 'double':
        # Double triangle
        X = np.asarray([[0,1,0],[-1,-1,0],[1,-1,0]], dtype=np.float64)
        X = np.concatenate([X, X*(-1,1,1)])
        tris = np.asarray([[0,1,2],[3,4,5]])
        segs = np.asarray([[[-2,y,y],[2,y,y]] for y in np.linspace(-1,1)], dtype=np.float64)
        clipped = clip_hidden(X[tris], segs)
        plot_segs(clipped)
