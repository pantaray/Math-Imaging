# myvec.py - Plotting tools for vector fields
# 
# Author: Stefan Fuertinger [stefan.fuertinger@gmx.at]
# Created: June 13 2012
# Last modified: <2017-09-21 12:41:45>

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

##########################################################################################
def myquiv(u,v):
    """
    Plots a 2D vector field using "sane" defaults for Matplotlib's `quiver`.

    Parameters
    ----------
    u : NumPy 2darray
        `x` components of the vector field `w(x,y) = (u(x,y),v(x,y))`. Note that `u` has to be 
        a 2D array of the same dimension as `v`. 
    v : NumPy 2darray
        `y` components of the vector field `w(x,y) = (u(x,y),v(x,y))`. Note that `v` has to be 
        a 2D array of the same dimension as `u`. 
       
    Returns
    -------
    Nothing : None 

    See also
    --------
    quiver : in the `Matplotlib Example Code Repository <http://matplotlib.org/examples/pylab_examples/quiver_demo.html>`_
    """

    # Check the input vector field
    checkfield(u,v,varnames = ["u","v"])

    # Now do something
    N  = u.shape[0]
    dN = min(N,16)
    plt.quiver(v[N-1:0:(-N/dN),0:N:(N/dN)],u[N-1:0:(-N/dN),0:N:(N/dN)],color="k")
    plt.axis("image")
    plt.axis("off")
    plt.draw()

    return

##########################################################################################
def makegrid(N,M=None,xmin=1,xmax=None,ymin=1,ymax=None):
    """
    Create an `M`-by-`N` grid on the 2D-domain `[xmin,xmax]`-by-`[ymin,ymax]`

    Parameters
    ----------
    N : int
        The number of grid-points in vertical (i.e. `y`-) direction
    M : int 
        The number of grid-points in horizontal (i.e. `x`-) direction. By default `M = N`
    xmin : float
        The left boundary of the (rectangular) domain. By default `xmin = 1` 
    xmax : float
        The right boundary of the (rectangular) domain. By default `xmax = N`
    ymin : float
        The lower boundary of the (rectangular) domain. By default `ymin = 1`
    ymax : float
        The upper boundary of the (rectangular) domain. By default `ymax = N`
    
    Returns
    -------
    x : NumPy 2darray
        2D grid array of `x`-values on the domain `[xmin,xmax]`-by-`[ymin,ymax]`
    y : NumPy 2darray
        2D grid array of `y`-values on the domain `[xmin,xmax]`-by-`[ymin,ymax]`

    Examples
    --------
    The call
 
    >>> x,y = makegrid(N)

    creates a square `[1,N]`-by-`[1,N]` grid given by

    >>> x
        array([[   1.,    1.,    1., ...,    1.],
               [   2.,    2.,    2., ...,    2.],
               ..., 
               [    N,     N,     N, ...,     N]])

    and

    >>> y
        array([[   1.,    2.,    3., ...,     N],
               [   1.,    2.,    3., ...,     N],
               ..., 
               [   1.,    2.,    3., ...,     N]])

    See also
    --------
    meshgrid : NumPy's `meshgrid <http://docs.scipy.org/doc/numpy/reference/generated/numpy.meshgrid.html>`_
    gengrid : found in the module `makeimg <makeimg.gengrid.html>`_
    """

    # Sanity checks
    scalarcheck(N,"N",kind=int,bounds=[1,np.inf])
    scalarcheck(M,"M",kind=int,bounds=[1,np.inf])
    scalarcheck(xmin,"xmin")
    scalarcheck(xmax,"xmax",bounds=[xmin,np.inf])
    scalarcheck(ymin,"ymin")
    scalarcheck(ymax,"ymax",bounds=[ymin,np.inf])

    # Compute stepsizes
    hx = (xmax - xmin)/(M-1)
    hy = (ymax - ymin)/(N-1)

    # Build 1D grid arrays
    x1 = xmin + hx*np.arange(0,M)
    y1 = ymin + hy*np.arange(0,N)

    # Build 2D grid arrays
    y,x = np.meshgrid(x1,y1)

    return x,y

##########################################################################################
def mygrid(u,v,x=None,y=None,rowstep=16,colstep=16,interpolation="lanczos"):
    """
    Plot a 2D vector field as deformed grid on a 2D lattice

    Parameters
    ----------
    u : NumPy 2darray
        `x` components of the vector field `w(x,y) = (u(x,y),v(x,y))`. Note that `u` has to be 
        a 2D array of the same dimension as `v`. 
    v : NumPy 2darray
        `y` components of the vector field `w(x,y) = (u(x,y),v(x,y))`. Note that `v` has to be 
        a 2D array of the same dimension as `u`. 
    x : NumPy 2darray
        2D grid array of `x`-values on the domain of `w(x,y)`. By default it is assumed that 
        `w` is defined on `[1,N]`-by-`[1,N]` 
    y : NumPy 2darray
        2D grid array of `y`-values on the domain of `w(x,y)`. By default it is assumed that 
        `w` is defined on `[1,N]`-by-`[1,N]` 
    rowstep : int
        Array row stride (step size) used to generate the grid. Default value is 16. 
    colstep : int
        Array column stride (step size) used to generate the grid. Default value is 16. 
    interpolation : str
        Interpolation to be used for plotting. Default value is "lanczos". 
        Recommended other values are "bilinear" or "nearest". See Matplotlib's 
        `imshow`-documentation for details. 
       
    Returns
    -------
    Nothing : None 

    See also
    --------
    imshow : in the `Matplotlib documentation <http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow>`_
    """

    # Check the input vector field
    checkfield(u,v,varnames=["u","v"])
    (M,N) = u.shape

    # Check the grid
    if (x is None and y != None) or (x != None and y == None):
        print "WARNING: Only x- or y-grid-data specified - switching to default domain `[1,N]`-by-`[1,N]`..."
        x,y = makegrid(N,M=M)
    elif x == None:
        x,y = makegrid(N,M=M)
    else:
        checkfield(x,y,varnames=["x","y"])
        checkgrid(u,x,y)

    # Sanity checks
    scalarcheck(rowstep,"rowstep",kind="int",bounds=[1,u.shape[0]])
    scalarcheck(colstep,"colstep",kind="int",bounds=[1,u.shape[1]])
    if not isinstance(interpolation,(str,unicode)):
        raise TypeError('Input `interpolation` has to be a string!')

    # Create lattice
    wires = np.ones(u.shape)
    wires[0:-1:rowstep,:] = 0
    wires[:,0:-1:colstep] = 0

    # Apply vector field
    zipyx  = zip(y.flatten(1),x.flatten(1))
    wiresr = griddata(zipyx,wires.flatten(1),(y+v,x+u),method="linear",fill_value=1)
    
    # Plot it
    plt.imshow(wiresr,cmap="gray",interpolation=interpolation)
    plt.axis("image")
    plt.axis("off")
    plt.draw()

    return

##########################################################################################
def mywire(u,v,x=None,y=None,rowstep=1,colstep=1):
    """
    Plot a 2D vector field as 3D wire-frame using Matplotlib's `plot_wireframe`

    Parameters
    ----------
    u : NumPy 2darray
        `x` components of the vector field `w(x,y) = (u(x,y),v(x,y))`. Note that `u` has to be 
        a 2D array of the same dimension as `v`. 
    v : NumPy 2darray
        `y` components of the vector field `w(x,y) = (u(x,y),v(x,y))`. Note that `v` has to be 
        a 2D array of the same dimension as `u`. 
    x : NumPy 2darray
        2D grid array of `x`-values on the domain of `w(x,y)`. By default it is assumed that 
        `w` is defined on `[1,N]`-by-`[1,N]` 
    y : NumPy 2darray
        2D grid array of `y`-values on the domain of `w(x,y)`. By default it is assumed that 
        `w` is defined on `[1,N]`-by-`[1,N]` 
    rowstep : int
        Array row stride (step size) used to generate the wire-frame plot (see `plot_wireframe`'s 
        documentation for details). Default value is 1. 
    colstep : int
        Array column stride (step size) used to generate the wire-frame plot (see `plot_wireframe`'s 
        documentation for details). Default value is 1. 
       
    Returns
    -------
    Nothing : None 

    See also
    --------
    plot_wireframe : in the `Matplotlib documentation <http://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D.plot_wireframe>`_
    """

    # Check the input vector field
    checkfield(u,v,varnames=["u","v"])
    (M,N) = u.shape

    # Check the grid
    if (x is None and y != None) or (x != None and y == None):
        print "WARNING: Only x- or y-grid-data specified - switching to default domain `[1,N]`-by-`[1,N]`..."
        x,y = makegrid(N,M=M)
    elif x == None:
        x,y = makegrid(N,M=M)
    else:
        checkfield(x,y,varnames=["x","y"])
        checkgrid(u,x,y)

    # Sanity checks
    scalarcheck(rowstep,"rowstep",kind="int",bounds=[1,u.shape[0]])
    scalarcheck(colstep,"colstep",kind="int",bounds=[1,u.shape[1]])

    # Draw the deformed grid
    ax = plt.gca(projection='3d')
    ax.plot_wireframe(x,y,u+v,color="black",rstride=rowstep,cstride=colstep)
    ax.set_axis_off()
    plt.draw()

    return

##########################################################################################
def checkfield(u,v,varnames=["u","v"]):
    """
    Perform sanity checks on the input vector field
    """
    
    for nk, arr in enumerate([u,v]):
        varname = varnames[nk]
        try:
            sha = arr.shape
        except:
            raise TypeError('Input `'+varname+'` must be a NumPy array, not '+type(arr).__name__+'!')
        if len(sha) != 2:
            raise ValueError('Input `'+varname+'` must be a `M`-by-`N` NumPy array')
        if (min(sha)==1) or (sha[0]!=sha[1]):
            raise ValueError('Input `'+varname+'` must be a `M`-by-`N` NumPy array!')
        if not plt.is_numlike(arr) or not np.isreal(arr).all():
            raise TypeError('Input `'+varname+'` must be a real-valued `M`-by-`N` NumPy array!')
        if np.isfinite(arr).min() == False:
            raise ValueError('Input `'+varname+'` must be a real valued NumPy array without Infs or NaNs!')
    if u.shape[0] != v.shape[0] or u.shape[1] != v.shape[1]:
        raise ValueError("Inputs `"+varnames[0]+"` and `"+varnames[1]+"` must have the same dimension!")
        
    return

##########################################################################################
def checkgrid(u,x,y):
    """
    Perform sanity checks on the grid
    """

    if x.shape[0] != u.shape[0] or x.shape[1] != u.shape[1]:
        raise ValueError("Grid and vector field must have the same dimensions!")
    if u.shape[0] != y.shape[0] or u.shape[1] != y.shape[1]:
        raise ValueError("Grid and vector field must have the same dimensions!")
        
    return

##########################################################################################
def scalarcheck(val,varname,kind=None,bounds=None):
    """
    Local helper function performing sanity checks on scalars
    """

    if not np.isscalar(val) or not plt.is_numlike(val) or not np.isreal(val).all():
        raise TypeError("Input `"+varname+"` must be a real scalar!")
    if not np.isfinite(val):
        raise TypeError("Input `"+varname+"` must be finite!")

    if kind == 'int':
        if (round(val) != val):
            raise ValueError("Input `"+varname+"` must be an integer!")

    if bounds is not None:
        if val < bounds[0] or val > bounds[1]:
            raise ValueError("Input scalar `"+varname+"` must be between "+str(bounds[0])+" and "+str(bounds[1])+"!")
