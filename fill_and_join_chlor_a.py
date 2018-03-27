#!/usr/bin/env python

import argparse
import netCDF4
from numba import jit
import numpy
import os
import pickle
import scipy.sparse
import scipy.sparse.linalg
import sys
import time

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  """

  # Arguments
  parser = argparse.ArgumentParser(description=
      """
      fill_and_join_chlor_a.py fills missing data in chlorophyll monthly climatology and joins into single file
      """,
      epilog='Written by A.Adcroft, 2018.')
  parser.add_argument('new_file', type=str,
      help="""Filename for new filled and joined dataset.""")
  parser.add_argument('bathy_file', type=str,
      help="""Filename for ocean bathymetry.""")
  parser.add_argument('files', type=str, nargs='+',
      help="""Filenames of sequential data.""")
  parser.add_argument('-bv','--bathy_var', type=str, default='depth',
      help="""Variable for depth in ocean bathymetry.""")
  parser.add_argument('-cv','--chlor_var', type=str, default='chlor_a',
      help="""Variable for chlorophyll in sequential files.""")
  parser.add_argument('-p','--progress', action='store_true',
      help="""Show progress.""")
  parser.add_argument('-q','--quiet', action='store_true',
      help="""Disable informational messages.""")

  return parser.parse_args()

def info(msg):
  """Prints an informational message with trailing ... and no newline"""
  print(msg + ' ...', end='')
  sys.stdout.flush()
  return time.time()

def end_info(tic):
  """Closes the informational line"""
  print(' done in %.3fs.'%(time.time() - tic) )

def main(args):
  """
  Does everything.
  """

  start_time = time.time()
  pickle_file = 'pickle.regrid_runoff_A'

  # Open bathymetry grid
  if args.progress: tic = info('Reading ocean bathymetry')
  ocn_depth = netCDF4.Dataset(args.bathy_file).variables[args.bathy_var]
  ocn_lat = netCDF4.Dataset(args.bathy_file).variables[ocn_depth.dimensions[0]][:]
  ocn_lon = netCDF4.Dataset(args.bathy_file).variables[ocn_depth.dimensions[1]][:]
  ocn_nj, ocn_ni = ocn_depth.shape
  ocn_depth = ocn_depth[:]
  if args.progress: end_info(tic)

  if not args.quiet: print('Bathymetry grid shape is %i x %i.'%(ocn_nj, ocn_ni))

  # Read river grid
  if args.progress: tic = info('Reading chlorophyll grid')
  nc = netCDF4.Dataset(args.files[0])
  chl = nc.variables[args.chlor_var]
  chl_nj, chl_ni = chl.shape
  chl_lon = nc.variables[chl.dimensions[1]][:]
  chl_lat = nc.variables[chl.dimensions[0]][::-1]
  if args.progress: end_info(tic)

  if not args.quiet: print('Chlorophyll grid shape is %i x %i.'%(chl_nj, chl_ni))

  # Create mask
  if args.progress: tic = info('Making up chlorophyll mask')
  yres, xres = ocn_nj / chl_nj, ocn_ni / chl_ni
  if xres != int(xres): raise Exception('i-size of bathymetry and chlorophyll data are not related.')
  if yres != int(yres): raise Exception('i-size of bathymetry and chlorophyll data are not related.')
  yres, xres = int(yres), int(xres)
  ocn_mask = numpy.zeros( ocn_depth.shape )
  ocn_mask[ ocn_depth[:,:]<0 ] = 1
  del ocn_depth
  ocn_mask = ocn_mask.reshape((chl_nj,yres,chl_ni,xres)).sum(axis=-1).sum(axis=1)
  ocn_mask[ ocn_mask>0 ] = 1
  ocn_lat = ocn_lat.reshape((chl_nj,yres)).mean(axis=1)
  ocn_lon = ocn_lon.reshape((chl_ni,xres)).mean(axis=1)
  if numpy.max( numpy.abs( ocn_lon - chl_lon ) ) > 1.e-4: raise Exception('Inconsistent longitudes!')
  if numpy.max( numpy.abs( ocn_lat - chl_lat ) ) > 1.e-4: raise Exception('Inconsistent latitudes!')
  if args.progress: end_info(tic)

  if not args.quiet: print('# of wet points = %i.'%(ocn_mask.sum()))

  # Remove lakes
  if args.progress: tic = info('Finding contiguous ocean cells')
  j,i = numpy.argmin( numpy.abs(ocn_lat+50) ), numpy.argmin( numpy.abs(ocn_lon-0) )
  ocn_mask = ice9it(int(j), int(i), ocn_mask)
  del j,i
  if args.progress: end_info(tic)

  if not args.quiet: print('# of wet points after ice-9 = %i.'%(ocn_mask.sum()))

  # Create file
  mo = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  mid, b = numpy.zeros(12), 0
  for n in range(12):
    mid[n] = b + mo[n]/2
    b += mo[n]
  del mo, b

  new_file = netCDF4.Dataset(args.new_file, 'w', 'clobber', format="NETCDF3_64BIT_OFFSET")
  # Global attributes
  for a in nc.ncattrs():
    new_file.setncattr(a, nc.getncattr(a))
  # Dimensions
  for d in nc.dimensions:
    if d not in ('rgb', 'eightbitcolor'):
      new_file.createDimension(d,len(nc.dimensions[d]))
  new_file.createDimension('time',None)
  for v in nc.variables:
    if v not in (args.chlor_var, 'palette'):
      nv = new_file.createVariable(v, nc.variables[v].dtype,  nc.variables[v].dimensions)
      for a in nc.variables[v].ncattrs():
        nv.setncattr(a,nc.variables[v].getncattr(a))
  tim = new_file.createVariable('time','f4',('time',))
  tim.modulo = ' '
  tim.units = 'days since 0001-01-01 00:00:00'
  tim.calendar = 'noleap'
  tim.long_name = 'Time'
  v = args.chlor_var
  nv = new_file.createVariable(v, nc.variables[v].dtype,  ('time','lat','lon',) )
  for a in nc.variables[v].ncattrs():
    nv.setncattr(a,nc.variables[v].getncattr(a))
  for v in nc.variables:
    if v not in (args.chlor_var, 'palette'):
      if v == 'lat':
        new_file.variables[v][:] = nc.variables[v][::-1]
      else:
        new_file.variables[v][:] = nc.variables[v][:]
  for n in range(len(args.files)):
    print(mid[n], args.files[n])
    chlor_a = netCDF4.Dataset(args.files[n]).variables[args.chlor_var][::-1,:]
    ndata = fill_missing_data(chlor_a, ocn_mask, chl_lat)
    nv[n] = ndata[:]
    tim[n] = mid[n]
  new_file.close()


def ice9it(j0,i0,depth):
  # Iterative implementation of "ice 9"
  wetMask = 0*depth
  (nj,ni) = wetMask.shape
  stack = set()
  stack.add( (j0,i0) )
  while stack:
    (j,i) = stack.pop()
    if (wetMask[j,i]>0) or (depth[j,i] == 0): continue
    wetMask[j,i] = 1
    if i>0: stack.add( (j,i-1) )
    else: stack.add( (j,ni-1) )
    if i<ni-1: stack.add( (j,i+1) )
    else: stack.add( (0,j) )
    if j>0: stack.add( (j-1,i) )
    if j<nj-1: stack.add( (j+1,i) )
    else: stack.add( (j,ni-1-i) )
  return wetMask

def fill_missing_data(idata, mask, lat, fast=True, verbose=True, maxiter=0, debug=False, stabilizer=1.e-14):
    """
    Returns data with masked values objectively interpolated except where mask==0.
    
    Arguments:
    data - numpy.ma.array with mask==True where there is missing data or land.
    mask - numpy.array of 0 or 1, 0 for land, 1 for ocean.
    
    Returns a numpy.ma.array.
    """
    nj,ni = idata.shape
    fdata = idata.filled(0.) # Working with an ndarray is faster than working with a masked array
    if debug:
        plt.figure(); plt.pcolormesh(mask); plt.title('mask'); plt.colorbar();
        plt.figure(); plt.pcolormesh(idata.mask); plt.title('idata.mask'); plt.colorbar();
        plt.figure(); plt.pcolormesh(idata); plt.title('idata'); plt.colorbar();
        plt.figure(); plt.pcolormesh(idata.filled(3.)); plt.title('idata.filled'); plt.colorbar();
        plt.figure(); plt.pcolormesh(idata.filled(3.)); plt.title('fdata'); plt.colorbar();
    missing_j, missing_i = numpy.where( idata.mask & (mask>0) )
    n_missing = missing_i.size
    if verbose:
        print('Data shape: %i x %i = %i with %i missing values'%(nj, ni, nj*ni, numpy.count_nonzero(idata.mask)))
        print('Mask shape: %i x %i = %i with %i land cells'%(mask.shape[0], mask.shape[1],
                                                                 numpy.prod(mask.shape), numpy.count_nonzero(1-mask)))
        print('Data has %i missing values in ocean'%(n_missing))
        print('Data range: %g .. %g '%(idata.min(),idata.max()))
    # ind contains column of matrix/row of vector corresponding to point [j,i]
    ind = numpy.zeros( fdata.shape, dtype=int ) - int(1e6)
    ind[missing_j,missing_i] = numpy.arange( n_missing )
    if verbose: print('Building matrix')
    A = scipy.sparse.lil_matrix( (n_missing, n_missing) )
    b = numpy.zeros( (n_missing) )
    ld = numpy.zeros( (n_missing) )
    if fast:
        if verbose: print('Constructing diagonals')
        # North
        jm,im = numpy.minimum(missing_j+1,nj-1),missing_i
        si = ind[jm,im]
        A[(si>=0) & (missing_j<nj-1),si[(si>=0) & (missing_j<nj-1)]] = 1.
        ll = (si<0) & (missing_j<nj-1)
        b[ll] -= fdata[jm, im][ll] * mask[jm,im][ll]
        ld[(missing_j<nj-1)] -= mask[jm, im][(missing_j<nj-1)]
        # South
        jm,im = numpy.maximum(missing_j-1,0),missing_i
        si = ind[jm,im]
        A[(si>=0) & (missing_j>0),si[(si>=0) & (missing_j>0)]] = 1.
        ll = (si<0) & (missing_j>0)
        b[ll] -= fdata[jm, im][ll] * mask[jm,im][ll]
        ld[(missing_j>0)] -= mask[jm, im][(missing_j>0)]
        # East
        jm,im = missing_j,numpy.mod(missing_i+1, ni)
        si = ind[jm,im]
        A[si>=0,si[si>=0]] = 1. /(numpy.cos(lat[jm[si>=0]]*numpy.pi/180)**2)
        ll = (si<0)
        b[ll] -= fdata[jm, im][ll] * mask[jm,im][ll] /(numpy.cos(lat[jm[si<0]]*numpy.pi/180)**2)
        ld[:] -= mask[jm, im] /(numpy.cos(lat[jm]*numpy.pi/180)**2)
        # West
        jm,im = missing_j,numpy.mod(missing_i+ni-1, ni)
        si = ind[jm,im]
        A[si>=0,si[si>=0]] = 1. /(numpy.cos(lat[jm[si>=0]]*numpy.pi/180)**2)
        ll = (si<0)
        b[ll] -= fdata[jm, im][ll] * mask[jm,im][ll]/(numpy.cos(lat[jm[si<0]]*numpy.pi/180)**2)
        ld[:] -= mask[jm, im] /(numpy.cos(lat[jm]*numpy.pi/180)**2)
    else:
        A[range(n_missing),range(n_missing)] = 0.
        if verbose: print('Looping over cells')
        for n in range(n_missing):
            j,i = missing_j[n],missing_i[n]
            im1 = ( i + ni - 1 ) % ni
            ip1 = ( i + 1 ) % ni
            jm1 = max( j-1, 0)
            jp1 = min( j+1, nj-1)
            if j>0 and mask[jm1,i]>0:
                ld[n] -= 1.
                ij = ind[jm1,i]
                if ij>=0:
                    A[n,ij] = 1.
                else:
                    b[n] -= fdata[jm1,i]
            if mask[j,im1]>0:
                ld[n] -= 1./(numpy.cos(lat[j]*numpy.pi/180)**2)
                ij = ind[j,im1]
                if ij>=0:
                    A[n,ij] = 1./(numpy.cos(lat[j]*numpy.pi/180)**2)
                else:
                    b[n] -= fdata[j,im1]/(numpy.cos(lat[j]*numpy.pi/180)**2)
            if mask[j,ip1]>0:
                ld[n] -= 1./(numpy.cos(lat[j]*numpy.pi/180)**2)
                ij = ind[j,ip1]
                if ij>=0:
                    A[n,ij] = 1./(numpy.cos(lat[j]*numpy.pi/180)**2)
                else:
                    b[n] -= fdata[j,ip1]/(numpy.cos(lat[j]*numpy.pi/180)**2)
            if j<nj-1 and mask[jp1,i]>0:
                ld[n] -= 1.
                ij = ind[jp1,i]
                if ij>=0:
                    A[n,ij] = 1.
                else:
                    b[n] -= fdata[jp1,i]
    if debug:
        tmp = numpy.zeros((nj,ni)); tmp[ missing_j, missing_i ] = b
        plt.figure(); plt.pcolormesh(tmp); plt.title('b (initial)'); plt.colorbar();
    # Set leading diagonal
    b[ld>=0] = 0.
    A[range(n_missing),range(n_missing)] = ld - stabilizer
    if debug:
        tmp = numpy.zeros((nj,ni)); tmp[ missing_j, missing_i ] = b
        plt.figure(); plt.pcolormesh(tmp); plt.title('b (final)'); plt.colorbar();
        tmp = numpy.ones((nj,ni)); tmp[ missing_j, missing_i ] = A.diagonal()
        plt.figure(); plt.pcolormesh(tmp); plt.title('A[i,i]'); plt.colorbar();
    if verbose: print('Matrix constructed')
    A = scipy.sparse.csr_matrix(A)
    if verbose: print('Matrix converted')
    new_data = numpy.ma.array( fdata, mask=(mask==0))
    if maxiter is None:
        x,info = scipy.sparse.linalg.bicg(A, b)
    elif maxiter==0:
        x = scipy.sparse.linalg.spsolve(A, b)
    else:
        x,info = scipy.sparse.linalg.bicg(A, b, maxiter=maxiter)
    if verbose: print('Matrix inverted')
    new_data[missing_j,missing_i] = x
    return new_data

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__':
  args = parseCommandLine()
  main(args)
