#!/usr/bin/env python
"""
Name: paraview_vtk_transects_extractor.py
Authors: Xylar Asay-Davis
Date: 04/17/2016

This script is used to extract a transect (a vertical slice through a 3D
field) from a time series of NetCDF files, using transects stored by
a NetCDF file produced using the MaskCreator tool.  The transects are
stored in VTK files for plotting in paraview.

It can extract a field across multiple files by passing in a regular expression
(or a semicolon-separated list of expressions) for the filename pattern. As an 
example, one can run the script using:

`./paraview_vtk_transect_extractor.py -v temperature -f "hist.comp.*.nc"`

The syntax is similar to paraview_vtk_field_extractor.py, except
that only fields on cells are currently supported and the following
additional flags are supported:
--transects_file: a NetCDF file produce by the MaskCreator that
  contains one or more transects

--transect_dim: the depth dimension for transects (e.g. nVertLevels
  in MPAS-Ocean)

--min_level_cell: a cell index array specifying the first valid
  index in each column (0 by default)

--max_level_cell: a cell index array specifiying the last valid
  index in each column (transect_dim-1 by default), e.g. maxLevelCell
  in MPAS-Ocean

--layer_thickness: the name of the field (optionally with a
  minus sign in front) where layer thicknesses are stored, used to
  compute the zInterface field that can be used in Paraview as
  the vertical coordinate for transects.

--z_min: the name of a field (optionally with a minus sign in front)
  specifying the depth of the lower interface (in index space) of the
  first valid layer on cells, used to compute zInterface

--z_max: the name of a field (optionally with a minus sign in front)
  specifying the depth of the upper interface (in index space) of the
  last valid layer on cells, used to compute zInterface

The script averages layerThickness from cell centers to edges (or used
the layerThickness of the valid cell for edges and layers where only
neighboring cell is valid, e.g. at stair steps in the topography).
Then, a cumulative sum of layer thickness is computed starting at either
z_min or z_max (z_min is used if both are specified) and stored in
zInterface for use in Paraview.


Requirements:
This script requires access to the following non standard modules:
pyevtk
netCDF4
numpy

Optional modules:
progressbar
"""
import os
import numpy as np

from netCDF4 import Dataset as NetCDFFile
import argparse

import utils

try:
    from progressbar import ProgressBar, Percentage, Bar, ETA
    use_progress_bar = True
except:
    use_progress_bar = False
    
class Transect(object):
    def __init__(self, cellIndices, mesh_file, nLevels, minLevelCell, maxLevelCell, dtype):

        self.cellIndices = cellIndices
        self.nLevels = nLevels

        nCells = len(cellIndices)
        self.nCells = nCells
        self.nCellsAndEdges = 2*nCells-1
        
        self.minLevelCell = minLevelCell
        self.maxLevelCell = maxLevelCell

        cellMask = np.ones((nCells,nLevels), bool)
        for iLevel in range(nLevels):
            if minLevelCell is not None:
                cellMask[:,iLevel] = np.logical_and(cellMask[:,iLevel], iLevel >= minLevelCell)
            if maxLevelCell is not None:
                cellMask[:,iLevel] = np.logical_and(cellMask[:,iLevel], iLevel <= maxLevelCell)

        edgeMask = np.logical_or(cellMask[0:-1,:], cellMask[1:,:])
        
        
        self.mask = np.zeros((self.nCellsAndEdges,self.nLevels), dtype=bool)
        self.mask[0::2,:] = cellMask
        self.mask[1::2,:] = edgeMask


        cellsOnCell = mesh_file.variables['cellsOnCell'][cellIndices,:]-1
        edgesOnCell = mesh_file.variables['edgesOnCell'][cellIndices,:]-1

        edgeIndices = np.zeros(nCells-1,int)
        for indexInTransect in range(nCells-1):
            iCellNext = cellIndices[indexInTransect+1]
            iNeighbor = np.nonzero(cellsOnCell[indexInTransect,:] == iCellNext)[0][0]
            edgeIndices[indexInTransect] = edgesOnCell[indexInTransect,iNeighbor]


        x = np.zeros(2*nCells-1, dtype=dtype)
        y = np.zeros(2*nCells-1, dtype=dtype)
        z = np.zeros(2*nCells-1, dtype=dtype)
        transectCellIndices = np.zeros(2*nCells-2, dtype=int)

        x[0::2] = mesh_file.variables['xCell'][cellIndices]
        x[1::2] = mesh_file.variables['xEdge'][edgeIndices]
        y[0::2] = mesh_file.variables['yCell'][cellIndices]
        y[1::2] = mesh_file.variables['yEdge'][edgeIndices]
        z[0::2] = mesh_file.variables['zCell'][cellIndices]
        z[1::2] = mesh_file.variables['zEdge'][edgeIndices]
        transectCellIndices[0::2] = np.arange(0,nCells-1)
        transectCellIndices[1::2] = np.arange(1,nCells)

        pointCellIndices = []
        pointLevelIndices = []
        pointHorizIndices = []
        pointVertIndices = []
        X = []
        Y = []
        Z = []
        for index in range(self.nCellsAndEdges-1):
            iCell = transectCellIndices[index]
            for iLevel in range(nLevels):
                if not (self.mask[index,iLevel] and self.mask[index+1,iLevel]):
                    continue
                horizIndices = [index, index+1, index+1, index]
                vertIndices = [iLevel, iLevel, iLevel+1, iLevel+1]
                X.extend([x[i] for i in horizIndices])
                Y.extend([y[i] for i in horizIndices])
                Z.extend([z[i] for i in horizIndices])
                pointCellIndices.extend([iCell,iCell,iCell,iCell])
                pointLevelIndices.extend([iLevel,iLevel,iLevel,iLevel])
                pointHorizIndices.extend(horizIndices)
                pointVertIndices.extend(vertIndices)

        self.nPoints = len(X)
        self.nPolygons = self.nPoints/4
        
        self.points = (np.array(X, dtype=dtype), np.array(Y, dtype=dtype), np.array(Z, dtype=dtype))
        self.pointCellIndices = np.array(pointCellIndices, dtype=int)
        self.pointLevelIndices = np.array(pointLevelIndices, dtype=int)
        self.pointHorizIndices = np.array(pointHorizIndices, dtype=int)
        self.pointVertIndices = np.array(pointVertIndices, dtype=int)
        
        self.connectivity = np.arange(self.nPoints)
        self.offsets = 4 + 4*np.arange(self.nPolygons)

    def computeZInterface(self, layerThicknessCell, zMinCell, zMaxCell, dtype):
        
        cellMask = self.mask[0::2,:]
        edgeMask = self.mask[1::2,:]
        
        # first average layerThickness to transect points (cells and edges both)
        layerThicknessEdge = (cellMask[0:-1,:]*layerThicknessCell[0:-1,:]+ cellMask[1:,:]*layerThicknessCell[1:,:])
        denom = (1.0*cellMask[0:-1,:] + 1.0*cellMask[1:,:])
        layerThicknessEdge[edgeMask] /= denom[edgeMask]

        layerThickness = np.zeros((self.nCellsAndEdges,self.nLevels), dtype=dtype)
        layerThickness[0::2,:] = layerThicknessCell
        layerThickness[1::2,:] = layerThicknessEdge

        zInterface = np.zeros((self.nCellsAndEdges,self.nLevels+1))
        for iLevel in range(self.nLevels):
            zInterface[:,iLevel+1] = (zInterface[:,iLevel]
                + self.mask[:,iLevel]*layerThickness[:,iLevel])

        # work your way from either the min or the max to compute zInterface
        if zMinCell is not None:
            minLevelCell = self.minLevelCell.copy()
            minLevelCell[minLevelCell < 0] = self.nLevels-1
            zOffsetCell = zMinCell - zInterface[np.arange(0,self.nCellsAndEdges,2),minLevelCell]
            levelEdge = np.maximum(minLevelCell[0:-1], minLevelCell[1:])
        else:
            maxLevelCell = self.maxLevelCell.copy()
            zOffsetCell = zMaxCell - zInterface[np.arange(0,self.nCellsAndEdges,2),maxLevelCell+1]
            levelEdge = np.minimum(maxLevelCell[0:-1], maxLevelCell[1:])

        for iLevel in range(self.nLevels+1):
            zInterface[0::2,iLevel] += zOffsetCell

        zOffsetEdge = (-zInterface[np.arange(1,self.nCellsAndEdges-1,2),levelEdge]
           + 0.5*(zInterface[np.arange(0,self.nCellsAndEdges-2,2),levelEdge]
           + zInterface[np.arange(2,self.nCellsAndEdges,2),levelEdge]))

        for iLevel in range(self.nLevels+1):
            zInterface[1::2,iLevel] += zOffsetEdge

        return zInterface[self.pointHorizIndices, self.pointVertIndices]


def build_transects( mesh_file, transects_file, time_series_file, 
                     transect_dim_name, min_level_name, max_level_name, 
                     output_32bit ):#{{{

    if output_32bit:
        outType = 'f4'
    else:
        outType = 'f8'

    nTransects = len(transects_file.dimensions['nTransects'])
    nTransectLevels = len(time_series_file.dimensions[transect_dim_name])

    if (mesh_file is not None):
        nc_file = mesh_file
    else:
        nc_file = time_series_file        

    # list of points starts out empty; will add points from each transect
    transects = []

    for iTransect in range(nTransects):
        # get the global indices of this transect
        cellIndices = transects_file.variables['transectCellGlobalIDs'][iTransect,:]-1
        cellIndices = cellIndices[cellIndices >= 0]
        if min_level_name is not None:
            var = utils.get_var(min_level_name, mesh_file, time_series_file)
            minLevelCell = var[cellIndices]-1
        else:
            minLevelCell = None
        if max_level_name is not None:
            var = utils.get_var(max_level_name, mesh_file, time_series_file)
            maxLevelCell = var[cellIndices]-1
        else:
            maxLevelCell = None
        
        transect = Transect(cellIndices, nc_file, nTransectLevels, minLevelCell, 
                            maxLevelCell, outType)
                            
        transects.append(transect)    
        
    return transects

#}}}

def build_transects_time_series( local_time_indices, file_names, mesh_file, all_dim_vals,
                                 transect_dim_name, transects, variable_names,
                                 layer_thickness_name, z_min_name, z_max_name,
                                 output_32bit, combine_output, 
                                 append ):#{{{
        
    def read_transect_field(transect, field_name):
        if field_name[0] == '-':
            field_name = field_name[1:]
            sign = -1
        else:
            sign = 1
            
        field_var = variables[field_name]
        field_ndims = len(field_var.dimensions)
        if transect_dim_name in field_var.dimensions:
            field = np.zeros((transect.nCells,transect.nLevels), dtype=outType)
            
            for iCellTransect in range(transect.nCells):
                iCell = transect.cellIndices[iCellTransect]

                dim_vals = []
                for iDim in range(field_ndims):
                    dim = field_var.dimensions[iDim]
                    if dim == 'Time':
                        dim_vals.append(local_time_indices[time_index])
                    elif dim == 'nCells':
                        dim_vals.append(iCell)
                    elif dim == transect_dim_name:
                        dim_vals.append(np.arange(transect.nLevels))
                    else:
                        dim_vals.append(all_dim_vals[field_name][0])


                if field_ndims == 1:
                    field[iCellTransect,:] = field_var[dim_vals[0]]
                elif field_ndims == 2:
                    field[iCellTransect,:]= field_var[dim_vals[0], dim_vals[1]]
                elif field_ndims == 3:
                    field[iCellTransect,:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2]]
                elif field_ndims == 4:
                    field[iCellTransect,:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3]]
                elif field_ndims == 5:
                    field[iCellTransect,:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3], dim_vals[4]]

        else:
            field = np.zeros((transect.nCells), dtype=outType)

            dim_vals = []
            for iDim in range(field_ndims):
                dim = field_var.dimensions[iDim]
                if dim == 'Time':
                    dim_vals.append(local_time_indices[time_index])
                elif dim == 'nCells':
                    dim_vals.append(transect.cellIndices)
                else:
                    print field_name, all_dim_vals
                    dim_vals.append(all_dim_vals[field_name][0])

            if field_ndims == 1:
                field[:] = field_var[dim_vals[0]]
            elif field_ndims == 2:
                field[:]= field_var[dim_vals[0], dim_vals[1]]
            elif field_ndims == 3:
                field[:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2]]
            elif field_ndims == 4:
                field[:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3]]
            elif field_ndims == 5:
                field[:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3], dim_vals[4]]


        return sign*field

    if len(variable_names) == 0:
        return

    if output_32bit:
        outType = 'float32'
    else:
        outType = 'float64'


    # Get dimension info to allocate the size of Colors
    time_series_file = NetCDFFile(file_names[0], 'r')

    nTransects = len(transects)
    nTimes = len(local_time_indices)

        
    nVars = len(variable_names)
    var_has_time_dim = np.zeros(nVars,bool)
    variables = {}
    for iVar in range(nVars):
        var_name = variable_names[iVar]
        if var_name == 'zInterface':
            # This is a special case that we're going to build ourselves
            variables[var_name] = None
            var_has_time_dim[iVar] = True
            continue
        
        if (mesh_file is not None) and (var_name in mesh_file.variables):
            var = mesh_file.variables[var_name]
            if 'Time' in var.dimensions:
                # we can't support time dependence in the mesh file
                var = time_series_file.variables[var_name]

        else:
            var = time_series_file.variables[var_name]

        variables[var_name] = var

        dims = var.dimensions

        var_has_time_dim[iVar] = 'Time' in dims
         
    time_series_file.close()
    
    any_var_has_time_dim = np.any(var_has_time_dim)

    try:
        os.makedirs('vtk_files')
    except OSError:
        pass

    if any_var_has_time_dim:
        try:
            os.makedirs('vtk_files/time_series')
        except OSError:
            pass
    else:
        # there is no point in combining output if no fields have Time dim
        combine_output = False
        nTimes = 1

    # Output time series
    if use_progress_bar:
        widgets = ['Writing time series: ', Percentage(), ' ', Bar(), ' ', ETA()]
        field_bar = ProgressBar(widgets=widgets, maxval=nTimes*nVars*nTransects).start()
    else:
        print "Writing time series...."
        
    pad = np.array(np.floor(np.log10(nTransects)),int)+1
    template = '%%0%dd'%pad

    for iTransect in range(nTransects):
        transect = transects[iTransect]
        
        if any_var_has_time_dim:
            if combine_output or np.all(var_has_time_dim):
                out_prefix = "transect_%s"%(template%iTransect)
            else:
                out_prefix = "timeDependentTransect_%s"%(template%iTransect)
            # start the pvd file
            pvd_file = utils.write_pvd_header('vtk_files', out_prefix)
            pvd_file.write('<Collection>\n')

        if not combine_output and not np.all(var_has_time_dim):
            out_prefix = "staticTransect_%s"%(template%iTransect)
            varIndices = np.arange(nVars)[var_has_time_dim == False]
            timeIndependentFile = utils.write_vtp_header('vtk_files', out_prefix, varIndices[0],
                                                         varIndices, variable_names,
                                                         all_dim_vals, transect.points,
                                                         transect.connectivity,
                                                         transect.offsets,
                                                         transect.nPoints,
                                                         transect.nPolygons,
                                                         outType,
                                                         cellData=False,
                                                         pointData=True)
        prev_file = ""
        for time_index in range(nTimes):
    
            if prev_file != file_names[time_index]:
                if prev_file != "":
                    time_series_file.close()
                time_series_file = NetCDFFile(file_names[time_index], 'r')
                prev_file = file_names[time_index]

            if any_var_has_time_dim:
                # write the header for the vtp file
                vtp_file_prefix = "time_series/%s.%d"%(out_prefix, time_index)
                file_name = 'vtk_files/%s.vtp'%(vtp_file_prefix)
                if append and os.path.exists(file_name):
                    continue
    
                if combine_output:
                    varIndices = np.arange(nVars)
                else:
                    varIndices = np.arange(nVars)[var_has_time_dim]
                timeDependentFile = utils.write_vtp_header('vtk_files', vtp_file_prefix, varIndices[0],
                                                           varIndices, variable_names,
                                                           all_dim_vals, transect.points,
                                                           transect.connectivity,
                                                           transect.offsets,
                                                           transect.nPoints,
                                                           transect.nPolygons,
                                                           outType,
                                                           cellData=False,
                                                           pointData=True)
    
                # add time step to pdv file
                pvd_file.write('<DataSet timestep="%d" group="" part="0"\n'%(time_index))
                pvd_file.write('\tfile="%s.vtp"/>\n'%(vtp_file_prefix))
            
            for iVar in range(nVars):
                var_name = variable_names[iVar]
                
                if var_name == 'zInterface':
                    # build a cum sum of layer thickness for transect vertical coord
                        
                    vtkFile = timeDependentFile
        
                    layerThickness = read_transect_field(transect, layer_thickness_name)
                    zMin = None
                    zMax = None
                    if z_min_name is not None:
                        zMin = read_transect_field(transect, z_min_name)
                    elif z_max_name is not None:
                        zMax = read_transect_field(transect, z_max_name)
                    else:
                        zMax = np.zeros(transect.nCells)
        
                    zInterface = transect.computeZInterface(layerThickness, zMin, zMax, dtype=outType)
                                    
                    vtkFile.appendData(zInterface)
                else:

                    field_var = variables[var_name]
    
                    has_time = 'Time' in field_var.dimensions
                    if not combine_output and not has_time and time_index > 0:
                        continue
                    
                    if has_time or combine_output:
                        vtkFile = timeDependentFile
                    else:
                        vtkFile = timeIndependentFile
                        
                    field = read_transect_field(transect, var_name)
                    
                    if(len(field.shape) == 2):
                        field = field[transect.pointCellIndices,transect.pointLevelIndices]
                    else:
                        field = field[transect.pointCellIndices]
    
                    vtkFile.appendData(field)
                    
                    if use_progress_bar:
                        field_bar.update((iTransect*nTimes + time_index)*nVars + iVar)
    
                    del field

                
            if any_var_has_time_dim:
                timeDependentFile.save()
                del timeDependentFile
    
            if time_index == 0 and not combine_output and not np.all(var_has_time_dim):
                timeIndependentFile.save()
                del timeIndependentFile
    
        time_series_file.close()

    
        if any_var_has_time_dim:
            # finish the pdv file
            pvd_file.write('</Collection>\n')
            pvd_file.write('</VTKFile>\n')

    if use_progress_bar:
        field_bar.finish()
#}}}


if __name__ == "__main__":
    if use_progress_bar:
        print " -- Using progress bars --"
    else:
        print " -- Progress bars are not available--"
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file_pattern", dest="filename_pattern", help="MPAS Filename pattern.", metavar="FILE", required=True)
    parser.add_argument("--transects_file", dest="transects_filename", help="MPAS transects filename.", required=True)
    parser.add_argument("--transect_dim", dest="transect_dim_name", help="A dimension for transect layers.", required=True)
    parser.add_argument("--layer_thickness", dest="layer_thickness_name", help="Variable for layer thickness, used to compute zInterface", required=True)
    parser.add_argument("--min_level_cell", dest="min_level_cell_name", help="Index array indicating the minimum valid layer in a cell (default is 0 for all cells)")
    parser.add_argument("--max_level_cell", dest="max_level_cell_name", help="Index array indicating the maximum valid layer in a cell (default is the transect_dim-1 for all cells)")
    parser.add_argument("--z_min", dest="z_min_name", help="Variable specifying the depth of the lower interface (in index space) of the first valid layer on cells, used to compute zInterface")
    parser.add_argument("--z_max", dest="z_max_name", help="Variable specifying the depth of the upper interface (in index space) of the last valid layer on cells, used to compute zInterface")
    parser.add_argument("-d", "--dim_list", dest="dimension_list", nargs="+", help="A list of dimensions and associated indices.")
    parser.add_argument("-m", "--mesh_file", dest="mesh_filename", help="MPAS Mesh filename. If not set, it will use the first file in the -f flag as the mesh file.")
    parser.add_argument("-v", "--variable_list", dest="variable_list", help="List of variables to extract ('all' for all variables, 'allOnCells' for all variables on cells, etc.)", metavar="VAR", required=True)
    parser.add_argument("-3", "--32bit", dest="output_32bit", help="If set, the vtk files will be written using 32bit floats.", action="store_true")
    parser.add_argument("-c", "--combine", dest="combine_output", help="If set, time-independent fields are written to each file along with time-dependent fields.", action="store_true")
    parser.add_argument("-a", "--append", dest="append", help="If set, only vtp files that do not already exist are written out.", action="store_true")
    args = parser.parse_args()

    if not args.output_32bit:
        use_32bit = False
    else:
        use_32bit = True

    (time_indices, time_file_names) = utils.setup_time_indices(args.filename_pattern)

    separate_mesh_file = True
    if not args.mesh_filename:
        args.mesh_filename = time_file_names[0]
        separate_mesh_file = False
        
    transects_file = NetCDFFile(args.transects_filename, 'r')

    # Setting dimension values:
    time_series_file = NetCDFFile(time_file_names[0], 'r')
    if separate_mesh_file:
        mesh_file = NetCDFFile(args.mesh_filename, 'r')
    else:
        mesh_file = None
    extra_dims = utils.parse_extra_dims(args.dimension_list, time_series_file, 
                                        mesh_file, max_index_count=1)
                                  
    (all_dim_vals, cellVars, vertexVars, edgeVars) = utils.setup_dimension_values_and_sort_vars(
            time_series_file, mesh_file, args.variable_list, 
            extra_dims,
            basic_dims=['nCells', 'Time', args.transect_dim_name],
            include_dims=['nCells'])
    
    if len(cellVars) == 0:
        print "No variables to extract."
        exit(0)
        
    cellVars.append('zInterface')
    all_dim_vals['zInterface'] = None
            
    print " -- Building transects --"
    
    transects = build_transects( mesh_file, transects_file, time_series_file,
                                 args.transect_dim_name, args.min_level_cell_name,
                                 args.max_level_cell_name, use_32bit )


    time_series_file.close()

    utils.summarize_extraction(args.mesh_filename, time_indices, 
                               cellVars, vertexVars, edgeVars, 
                               transects_file=args.transects_filename)

    print " -- Extracting cell fields on transects --"

    build_transects_time_series( time_indices, time_file_names, mesh_file,
                                 all_dim_vals, args.transect_dim_name, transects,
                                 cellVars, args.layer_thickness_name, 
                                 args.z_min_name, args.z_max_name, use_32bit,
                                 args.combine_output, args.append )
    if separate_mesh_file:
        mesh_file.close()

    transects_file.close()

# vim: set expandtab:
