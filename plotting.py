"""
Plotting routines for tracking.

Functions include:
* background
* hist
* tracks
"""

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib.mlab import *
from matplotlib.pyplot import *
import inout
import os
import pdb
import matplotlib.ticker as ticker
import op
import netCDF4 as netCDF
import tools

def background(grid=None):
    """
    Plot basic TXLA shelf background: coastline, bathymetry, meridians, etc
    Can optionally input grid (so it doesn't have to be loaded again)
    """

    matplotlib.rcParams.update({'font.size': 20})#,'font.weight': 'bold'})

    if grid is None:
        loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
        grid = inout.readgrid(loc)

    # Do plot   
    grid['basemap'].drawcoastlines()
    grid['basemap'].fillcontinents('0.8')
    grid['basemap'].drawparallels(np.arange(18, 35), dashes=(1, 0), linewidth=0.15, labels=[1, 0, 0, 0])
    grid['basemap'].drawmeridians(np.arange(-100, -80), dashes=(1, 0), linewidth=0.15, labels=[0, 0, 0, 1])
    hold('on')
    contour(grid['xr'], grid['yr'], grid['h'], np.hstack(([10,20],np.arange(50,500,50))), colors='lightgrey', linewidths=0.25)

    # Outline numerical domain
    plot(grid['xr'][0,:], grid['yr'][0,:], 'k:')
    plot(grid['xr'][-1,:], grid['yr'][-1,:], 'k:')
    plot(grid['xr'][:,0], grid['yr'][:,0], 'k:')
    plot(grid['xr'][:,-1], grid['yr'][:,-1], 'k:')


def hist(lonp, latp, fname, tind='final', which='contour', \
            bins=(40,40), N=10, grid=None, xlims=None, ylims=None):
    """
    Plot histogram of given track data at time index tind.

    Inputs:
        lonp,latp   Drifter track positions in lon/lat [time x ndrifters]
        fname       Plot name to save
        tind        (optional) Default is 'final', in which case the final
                    position of each drifter in the array is found
                    and plotted. Alternatively, a time index 
                    can be input and drifters at that time will be plotted.
                    Note that once drifters hit the outer numerical boundary,
                    they are nan'ed out so this may miss some drifters.
        which       (optional) 'contour', 'pcolor', 'hexbin', 'hist2d' 
                    for type of plot used. Default 'hexbin'.
        bins        (optional) Number of bins used in histogram. Default (15,25).
        N           (optional) Number of contours to make. Default 10.
        grid        (optional) grid as read in by inout.readgrid()
        xlims       (optional) value limits on the x axis
        ylims       (optional) value limits on the y axis

    Note: Currently assuming we are plotting the final location 
    of each drifter regardless of tind.
    """

    if grid is None:
        loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
        grid = inout.readgrid(loc)

    # Change positions from lon/lat to x/y
    xp, yp = grid['basemap'](lonp, latp)
    # Need to retain nan's since basemap changes them to values
    ind = np.isnan(lonp)
    xp[ind] = np.nan
    yp[ind] = np.nan

    fig = figure(figsize=(12,10))
    background(grid) # Plot coastline and such

    # pdb.set_trace()

    if tind == 'final':
        # Find final positions of drifters
        xpc, ypc = tools.find_final(xp, yp)
    elif is_numlike(tind):
        xpc = xp[:,tind]
        ypc = yp[:,tind]
    else: # just plot what is input if some other string
        xpc = xp.flatten()
        ypc = yp.flatten()

    if which == 'contour':

        # Info for 2d histogram
        H, xedges, yedges = np.histogram2d(xpc, ypc,
                                range=[[grid['xr'].min(), \
                                grid['xr'].max()], \
                                [grid['yr'].min(), \
                                grid['yr'].max()]],
                                bins=bins)
        # Contour Plot
        XE, YE = np.meshgrid(op.resize(xedges,0), op.resize(yedges,0))
        d = (H/H.sum())*100
        # # from http://matplotlib.1069221.n5.nabble.com/question-about-contours-and-clim-td21111.html
        # locator = ticker.MaxNLocator(50) # if you want no more than 10 contours
        # locator.create_dummy_axis()
        # locator.set_bounds(0,1)#d.min(),d.max())
        # levs = locator()
        con = contourf(XE, YE, d.T, N)#,levels=levs)#(0,15,30,45,60,75,90,105,120))
        con.set_cmap('YlOrRd')
        # Horizontal colorbar below plot
        cax = fig.add_axes([0.3725, 0.25, 0.48, 0.02]) #colorbar axes
        cb = colorbar(con, cax=cax, orientation='horizontal')
        cb.set_label('Final drifter location (percent)')

        # Save figure into a local directory called figures. Make directory if it doesn't exist.
        if not os.path.exists('figures'):
            os.makedirs('figures')

        savefig('figures/' + fname + 'histcon.png',bbox_inches='tight')
        # savefig('figures/' + fname + 'histcon.pdf',bbox_inches='tight')

    elif which == 'pcolor':

        # Info for 2d histogram
        H, xedges, yedges = np.histogram2d(xpc, ypc,
                                range=[[grid['xr'].min(), \
                                grid['xr'].max()], \
                                [grid['yr'].min(), \
                                grid['yr'].max()]],
                                bins=bins)
        # Pcolor plot
        p = pcolor(xedges, yedges, (H.T/H.sum())*100, cmap='YlOrRd')

        # Set x and y limits
        # pdb.set_trace()
        if xlims is not None:
            xlim(xlims)
        if ylims is not None:
            ylim(ylims)

        # Horizontal colorbar below plot
        cax = fig.add_axes([0.3775, 0.25, 0.48, 0.02]) #colorbar axes
        cb = colorbar(p, cax=cax, orientation='horizontal')
        cb.set_label('Final drifter location (percent)')

        # Save figure into a local directory called figures. Make directory if it doesn't exist.
        if not os.path.exists('figures'):
            os.makedirs('figures')

        savefig('figures/' + fname + 'histpcolor.png', bbox_inches='tight')
        # savefig('figures/' + fname + 'histpcolor.pdf',bbox_inches='tight')

    elif which == 'hexbin':

        # C with the reduce_C_function as sum is what makes it a percent
        C = np.ones(len(xpc))*(1./len(xpc))*100
        hb = hexbin(xpc, ypc, C=C, cmap='YlOrRd', gridsize=bins[0], 
                extent=(grid['xr'].min(), grid['xr'].max(), 
                grid['yr'].min(), grid['yr'].max()), 
                reduce_C_function=sum)

        # Set x and y limits
        # pdb.set_trace()
        if xlims is not None:
            xlim(xlims)
        if ylims is not None:
            ylim(ylims)

        # Horizontal colorbar below plot
        cax = fig.add_axes([0.3775, 0.25, 0.48, 0.02]) #colorbar axes
        cb = colorbar(cax=cax, orientation='horizontal')
        cb.set_label('Final drifter location (percent)')

        # pdb.set_trace()
        # Save figure into a local directory called figures. Make directory if it doesn't exist.
        if not os.path.exists('figures'):
            os.makedirs('figures')

        savefig('figures/' + fname + 'histhexbin.png', bbox_inches='tight')
        # savefig('figures/' + fname + 'histhexbin.pdf',bbox_inches='tight')

    elif which == 'hist2d':
        # pdb.set_trace()

        hist2d(xpc, ypc, bins=40, 
                range=[[grid['xr'].min(), grid['xr'].max()], 
                [grid['yr'].min(), grid['yr'].max()]], normed=True)
        set_cmap('YlOrRd')
        # Set x and y limits
        # pdb.set_trace()
        if xlims is not None:
            xlim(xlims)
        if ylims is not None:
            ylim(ylims)

        # Horizontal colorbar below plot
        cax = fig.add_axes([0.3775, 0.25, 0.48, 0.02]) #colorbar axes
        cb = colorbar(cax=cax,orientation='horizontal')
        cb.set_label('Final drifter location (percent)')

        # Save figure into a local directory called figures. Make directory if it doesn't exist.
        if not os.path.exists('figures'):
            os.makedirs('figures')

        savefig('figures/' + fname + 'hist2d.png',bbox_inches='tight')
        # savefig('figures/' + fname + 'histpcolor.pdf',bbox_inches='tight')


def tracks(lonp,latp,fname,grid=None):
    """
    Plot tracks as lines with starting points in green and ending points in red.

    Inputs:
        lonp,latp   Drifter track positions [time x ndrifters]
        fname       Plot name to save
    """
    # pdb.set_trace()
    if grid is None:
        loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
        grid = inout.readgrid(loc)

    # Change positions from lon/lat to x/y
    xp,yp = grid['basemap'](lonp,latp)
    # Need to retain nan's since basemap changes them to values
    ind = np.isnan(lonp)
    xp[ind] = np.nan
    yp[ind] = np.nan

    figure(figsize=(12,10))
    background(grid) # Plot coastline and such

    # pdb.set_trace()

    # Starting marker
    plot(xp[:,0],yp[:,0],'o',color='g',markersize=3,label='_nolegend_',alpha=0.4)

    # Plot tracks
    plot(xp.T,yp.T,'-',color='grey',linewidth=.2)

    # Find final positions of drifters
    xpc,ypc = tools.find_final(xp,yp)
    plot(xpc,ypc,'o',color='r',label='_nolegend_')
    # pdb.set_trace()

    # Legend, of sorts
    ax = gca()
    xtext = 0.45; ytext = 0.18;
    text(xtext, ytext, 'starting location', fontsize=16, color='green', 
        alpha=.8, transform = ax.transAxes)
    text(xtext, ytext-.03, 'track', fontsize=16, color='grey', transform = ax.transAxes)
    text(xtext, ytext-.03*2, 'ending location', fontsize=16, color='red', 
        transform = ax.transAxes)
    # xtext, ytext = grid['basemap'](-94,24) # text location
    # text(xtext,ytext,'starting location',fontsize=16,color='green',alpha=.8)
    # text(xtext,ytext-30000,'track',fontsize=16,color='grey')#,alpha=.8)
    # text(xtext,ytext-60000,'ending location',fontsize=16,color='red')#,alpha=.8)

    # # get psi mask from rho mask
    # # maskp = grid['mask'][1:,1:]*grid['mask'][:-1,1:]* \
 # #               grid['mask'][1:,:-1]*grid['mask'][:-1,:-1] 
    # # ind = maskp
    # # ind[ind==0] = np.nan
    # # plot(grid['xpsi']*ind,grid['ypsi']*ind,'k', \
    # #         (grid['xpsi']*ind).T,(grid['ypsi']*ind).T,'k')
    # plot(grid['xpsi'],grid['ypsi'],'k', \
    #       (grid['xpsi']).T,(grid['ypsi']).T,'k')

    # 16 is (lower) one that is near islands, 41 is higher one

    # show()

    # Save figure into a local directory called figures. Make directory if it doesn't exist.
    if not os.path.exists('figures'):
        os.makedirs('figures')

    savefig('figures/' + fname + 'tracks.png',bbox_inches='tight')
    # savefig('figures/' + fname + 'tracks.pdf',bbox_inches='tight')

def stream(lonp, latp, fname, which='contour', \
            grid=None, xlims=None, ylims=None):
    """
    UPDATE THIS
    Plot histogram of given track data at time index tind.

    Inputs:
        lonp,latp   Drifter track positions in lon/lat [time x ndrifters]
        fname       Plot name to save
        tind        (optional) Default is 'final', in which case the final
                    position of each drifter in the array is found
                    and plotted. Alternatively, a time index 
                    can be input and drifters at that time will be plotted.
                    Note that once drifters hit the outer numerical boundary,
                    they are nan'ed out so this may miss some drifters.
        which       (optional) 'contour', 'pcolor', 'hexbin', 'hist2d' 
                    for type of plot used. Default 'hexbin'.
        bins        (optional) Number of bins used in histogram. Default (15,25).
        N           (optional) Number of contours to make. Default 10.
        grid        (optional) grid as read in by inout.readgrid()
        xlims       (optional) value limits on the x axis
        ylims       (optional) value limits on the y axis

    Note: Currently assuming we are plotting the final location 
    of each drifter regardless of tind.
    """

    if grid is None:
        loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
        grid = inout.readgrid(loc)

    # Change positions from lon/lat to x/y
    xp, yp = grid['basemap'](lonp, latp)
    # Need to retain nan's since basemap changes them to values
    ind = np.isnan(lonp)
    xp[ind] = np.nan
    yp[ind] = np.nan

    fig = figure(figsize=(12,10))
    background(grid) # Plot coastline and such

    # pdb.set_trace()


    # C with the reduce_C_function as sum is what makes it a percent
    C = np.ones(len(xpc))*(1./len(xpc))*100
    hb = hexbin(xpc, ypc, C=C, cmap='YlOrRd', gridsize=40, 
            extent=(grid['xr'].min(), grid['xr'].max(), 
            grid['yr'].min(), grid['yr'].max()), 
            reduce_C_function=sum)

    # Set x and y limits
    # pdb.set_trace()
    if xlims is not None:
        xlim(xlims)
    if ylims is not None:
        ylim(ylims)

    # Horizontal colorbar below plot
    cax = fig.add_axes([0.3775, 0.25, 0.48, 0.02]) #colorbar axes
    cb = colorbar(cax=cax, orientation='horizontal')
    cb.set_label('Final drifter location (percent)')

    # pdb.set_trace()
    # Save figure into a local directory called figures. Make directory if it doesn't exist.
    if not os.path.exists('figures'):
        os.makedirs('figures')

    savefig('figures/' + fname + 'stream.png', bbox_inches='tight')
    # savefig('figures/' + fname + 'histhexbin.pdf',bbox_inches='tight')
