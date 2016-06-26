#!/usr/bin/env python
"""
Name: paraview_vtk_field_at_depth_extractor.py
Authors: Xylar Asay-Davis
Date: 06/26/2016

This script extracts fields at a specified depth (either given by
a constant depth value or a vertical index field on cells) as well
as the vertical location of the top and bottom interfaces of the
extracted layer.  

Currently, only fields on cells are supported. In paraview, fields
on cells will have constant values but will be interpolated over
rectangular steps between cells at edges.

Usage is similar to paraview_vtk_field_extractor.py. Additional flags are:
--vert_dim: specifies the vertical dimension of the model (e.g.
  nVertLevels in MPAS-Ocean)

--field_level_cell: specifies either a constant index value or the name 
  of an index array on cells that corresponds to the location of the 
  topography (e.g. maxLevelCell in MPAS-Ocean). If this flag is not specified, 
  the --depth flag should be used instead.

--depth: the constant depth at which fields should be extracted (if 
  field_level_cell is not specified)

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

Results are stored in topoFieldsOnCells or timeDependentTopoFieldsOnCells
and staticTopoFieldsOnCells (depending on the --combine flag), so as
not to conflict with fieldsOnCells, etc. that may be extracted separatedly.
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


def build_point_and_polygon_lists( nc_file, output_32bit, path ):#{{{

    cacheFileName = '%s/cache.nc'%path
    if os.path.exists(cacheFileName):
        # read in the points and connectivity info instead
        inCache = NetCDFFile(cacheFileName,'r')
        X = inCache.variables['X'][:]
        Y = inCache.variables['Y'][:]
        Z = inCache.variables['Z'][:]
        pointToCellMap = inCache.variables['pointToCellMap'][:]
        connectivity = inCache.variables['connectivity'][:]
        offsets = inCache.variables['offsets'][:]
        inCache.close()
        
        return ((X,Y,Z), connectivity, offsets, pointToCellMap)

    nCells = len(nc_file.dimensions['nCells'])
    nEdges = len(nc_file.dimensions['nEdges'])

    nEdgesOnCell = nc_file.variables['nEdgesOnCell'][:]
    verticesOnCell = nc_file.variables['verticesOnCell'][:,:]-1
    edgesOnCell = nc_file.variables['edgesOnCell'][:,:]-1
    xVertVar = nc_file.variables['xVertex']
    yVertVar = nc_file.variables['yVertex']
    zVertVar = nc_file.variables['zVertex']
    xCellVar = nc_file.variables['xCell']
    yCellVar = nc_file.variables['yCell']
    zCellVar = nc_file.variables['zCell']

    totalEdgesOnCells = np.sum(nEdgesOnCell)
    # cell center and 2 at each vertex
    nPoints = nCells+2*totalEdgesOnCells
    # one rectangle at each edge, a triangle between each edge and each cell
    # center and a triangle between the two points at each vertex and the
    # cell center
    nPolygons = 2*totalEdgesOnCells+nEdges

    if output_32bit:
        dtype = 'f4'
    else:
        dtype = 'f8'

    X = np.zeros(nPoints,dtype)
    Y = np.zeros(nPoints,dtype)
    Z = np.zeros(nPoints,dtype)
    pointToCellMap = np.zeros(nPoints,int)
    
    edgeVertToPointMap = -1*np.ones((nEdges,4),int)
    pointCountOnEdge = np.zeros((nEdges),int)
    
    offsets = np.zeros(nPolygons, dtype=int)
    # a trinangle per vertex and per edge on each cell and a rectangle per edge
    connectivity = np.zeros(6*totalEdgesOnCells+nEdges*4, dtype=int)
    
    
    if use_progress_bar:
        widgets = ['Build cell connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        cell_bar = ProgressBar(widgets=widgets, maxval=nCells).start()
    else:
        print "Build cell connectivity..."

    pointOffset = 0
    polyOffset = 0
    offsetIndex = 0
    
    for iCell in range(nCells):
        nVerts = nEdgesOnCell[iCell]
        iVerts = verticesOnCell[iCell,0:nVerts]
        iEdges = edgesOnCell[iCell,0:nVerts]
        
        # add point at the center of the cell
        X[pointOffset] = xCellVar[iCell]
        Y[pointOffset] = yCellVar[iCell]
        Z[pointOffset] = zCellVar[iCell]
        pointToCellMap[pointOffset] = iCell
        
        for offset in range(2):
            outIndices = pointOffset+1+offset+2*np.arange(nVerts)
            X[outIndices] = xVertVar[iVerts]
            Y[outIndices] = yVertVar[iVerts]
            Z[outIndices] = zVertVar[iVerts]
            pointToCellMap[outIndices] = iCell

        outIndices = pointOffset + 1 + np.mod(2*np.arange(nVerts)-1, 2*nVerts)
        edgeVertToPointMap[iEdges, pointCountOnEdge[iEdges]] = outIndices
        outIndices = pointOffset + 1 + 2*np.arange(nVerts)
        edgeVertToPointMap[iEdges, pointCountOnEdge[iEdges]+1] = outIndices
        pointCountOnEdge[iEdges] += 2
        
        # first vertex in every triangle is the cell center at index pointOffset
        outIndices = polyOffset + 3*np.arange(2*nVerts)
        connectivity[outIndices] = pointOffset
        # remaining vertices in each triangle are the remining points in order
        connectivity[outIndices+1] = pointOffset+1+np.arange(2*nVerts)
        connectivity[outIndices+2] = pointOffset+1+np.mod(1+np.arange(2*nVerts), 2*nVerts)
        # offset to the start of the next polygon
        offsets[offsetIndex:offsetIndex+2*nVerts] = polyOffset+3+3*np.arange(2*nVerts)

        
        pointOffset += 1 + 2*nVerts
        polyOffset += 6*nVerts
        offsetIndex += 2*nVerts
        
        if use_progress_bar:
            cell_bar.update(iCell)

    if use_progress_bar:
        cell_bar.finish()


    if use_progress_bar:
        widgets = ['Build edge connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        edge_bar = ProgressBar(widgets=widgets, maxval=nEdges).start()
    else:
        print "Build cell connectivity..."

    for iEdge in range(nEdges):
        connectivity[polyOffset] = edgeVertToPointMap[iEdge,0]
        connectivity[polyOffset+1] = edgeVertToPointMap[iEdge,1]
        if(pointCountOnEdge[iEdge] > 2):
            connectivity[polyOffset+2] = edgeVertToPointMap[iEdge,2]
            connectivity[polyOffset+3] = edgeVertToPointMap[iEdge,3]
        else:
            connectivity[polyOffset+2] = edgeVertToPointMap[iEdge,1]
            connectivity[polyOffset+3] = edgeVertToPointMap[iEdge,0]
        
        polyOffset += 4
        offsets[2*totalEdgesOnCells + iEdge] = polyOffset

        if use_progress_bar:
            edge_bar.update(iEdge)

    if use_progress_bar:
        edge_bar.finish()

    outCache = NetCDFFile(cacheFileName,'w')
    outCache.createDimension('nPoints',nPoints)
    outCache.createDimension('nPolygons',nPolygons)
    outCache.createDimension('nPolyVerts',len(connectivity))
    var = outCache.createVariable('X','float64',('nPoints'))
    var[:] = X
    var = outCache.createVariable('Y','float64',('nPoints'))
    var[:] = Y
    var = outCache.createVariable('Z','float64',('nPoints'))
    var[:] = Z
    
    var = outCache.createVariable('pointToCellMap','int32',('nPoints'))
    var[:] = pointToCellMap
    
    var = outCache.createVariable('connectivity','int32',('nPolyVerts'))
    var[:] = connectivity
    
    var = outCache.createVariable('offsets','int32',('nPolygons'))
    var[:] = offsets
    outCache.close()


    return ((X,Y,Z), connectivity, offsets, pointToCellMap)

#}}}


def build_field_time_series( local_time_indices, file_names, mesh_file, blocking, all_dim_vals,
                             variable_list, vertices, connectivity, offsets,
                             vert_dim_name, pointToCellMap, 
                             layer_thickness_name,
                             min_level_name, max_level_name,
                             z_min_name, z_max_name, 
                             field_level_name, depth,
                             output_32bit, combine_output, append, out_path ):#{{{

    def getVertInterpInfo(time_index):
        if min_level_name is not None:
            var = utils.get_var(min_level_name, mesh_file, time_series_file)
            minLevelCell = var[:]-1
        else:
            minLevelCell = None
    
        if max_level_name is not None:
            var = utils.get_var(max_level_name, mesh_file, time_series_file)
            maxLevelCell = var[:]-1
        else:
            maxLevelCell = None
            
        cellsOnEdge = utils.get_var('cellsOnEdge', mesh_file, time_series_file)[:,:]-1
        nEdgesOnCell = utils.get_var('nEdgesOnCell', mesh_file, time_series_file)[:]
        edgesOnCell = utils.get_var('edgesOnCell', mesh_file, time_series_file)[:,:]-1
            
        cellIndices = np.arange(nCells)
        (var_name, layerSign) = utils.get_field_sign(layer_thickness_name)
        var = utils.get_var(var_name, mesh_file, time_series_file)
        layerThickness = utils.read_field(field_var=var,
                                          extra_dim_vals=None,
                                          time_index=local_time_indices[time_index], 
                                          cell_indices=cellIndices,
                                          outType=outType,
                                          sign=layerSign,
                                          vert_dim_name=vert_dim_name,
                                          nLevels=nLevels)
        zMin = None
        zMax = None
        if z_min_name is not None:
            (var_name, sign) = utils.get_field_sign(z_min_name)
            var = utils.get_var(var_name, mesh_file, time_series_file)
            zMin = utils.read_field(field_var=var,
                                    extra_dim_vals=None, 
                                    time_index=local_time_indices[time_index], 
                                    cell_indices=cellIndices,
                                    outType=outType,
                                    sign=sign)
        elif z_max_name is not None:
            (var_name, sign) = utils.get_field_sign(z_max_name)
            var = utils.get_var(var_name, mesh_file, time_series_file)
            zMax = utils.read_field(field_var=var,
                                    extra_dim_vals=None, 
                                    time_index=local_time_indices[time_index], 
                                    cell_indices=cellIndices,
                                    outType=outType,
                                    sign=sign)
        else:
            zMax = np.zeros(nCells)

        (zInterfaceCell, zInterfaceEdge)  = utils.compute_zInterface(minLevelCell,
            maxLevelCell, layerThickness, zMin, zMax, dtype=outType,
            cellsOnEdge=cellsOnEdge)

        if depth is not None:
            valid = np.zeros(nCells,bool)
            zIndices = -1*np.ones(nCells,int)
            for zIndex in range(nLevels):
                if(layerSign == 1):
                    mask = np.logical_and(depth >= zInterfaceCell[:,zIndex], 
                                          depth < zInterfaceCell[:,zIndex+1])
                else:
                    mask = np.logical_and(depth <= zInterfaceCell[:,zIndex], 
                                          depth > zInterfaceCell[:,zIndex+1])
                valid[mask] = True
                zIndices[mask] = zIndex
        else:
            try:
                # first see if the field is an integer
                zIndices = int(field_level_name)*np.ones(nCells,int)
            except ValueError:
                # if not, it's the name of a variable
                var = utils.get_var(field_level_name, mesh_file, time_series_file)
                zIndices = var[:]-1
            valid = zIndices >= 0
        
        # extract the zInterfaces at the top and bottom of the layer at
        # cells and edges
        zTop = np.zeros(nPoints,float)
        zBot = np.zeros(nPoints,float)
        pointOffset = 0
        for iCell in range(nCells):
            if not valid[iCell]:
                continue
            zTop[pointOffset] = zInterfaceCell[iCell,zIndices[iCell]]
            zBot[pointOffset] = zInterfaceCell[iCell,zIndices[iCell]+1]
            
            nEdges = nEdgesOnCell[iCell]
            iEdges = edgesOnCell[iCell, 0:nEdges]
            outIndices = pointOffset + 1 + np.mod(2*np.arange(nEdges)-1,2*nEdges)
            zTop[outIndices] = zInterfaceEdge[iEdges,zIndices[iCell]]
            zBot[outIndices] = zInterfaceEdge[iEdges,zIndices[iCell]+1]
            outIndices = pointOffset + 1 + 2*np.arange(nEdges)
            zTop[outIndices] = zInterfaceEdge[iEdges,zIndices[iCell]]
            zBot[outIndices] = zInterfaceEdge[iEdges,zIndices[iCell]+1]
            pointOffset += 1 + 2*nEdges
            
        return (valid, zIndices, zTop, zBot)
        

    if len(variable_list) == 0:
        return

    if output_32bit:
        outType = 'float32'
    else:
        outType = 'float64'


    # Get dimension info to allocate the size of Colors
    time_series_file = NetCDFFile(file_names[0], 'r')

    if mesh_file is not None:
        # nGeomDim may not exist in time series file
        nCells = len(mesh_file.dimensions['nCells'])
    else:
        nCells = len(time_series_file.dimensions['nCells'])

    # Pre-compute the number of blocks
    nBlocks = 1 + nCells / blocking
        
    nPolygons = len(offsets)
    nPoints = len(vertices[0])
    nTimes = len(local_time_indices)
    nVars = len(variable_list)

    if vert_dim_name is not None:
        if vert_dim_name in time_series_file.dimensions:
            nLevels = len(time_series_file.dimensions[vert_dim_name])
        else:
            nLevels = len(mesh_file.dimensions[vert_dim_name])
    else:
        nLevels = None

    var_has_time_dim = np.zeros(nVars,bool)
    nHyperSlabs = 0
    for iVar in range(nVars):
        var_name = variable_list[iVar]
        if var_name in ['zInterfaceTop','zInterfaceBot']:
            # This is a special case that we're going to build ourselves
            var_has_time_dim[iVar] = True
            continue
            

        if (mesh_file is not None) and (var_name in mesh_file.variables):
            var = mesh_file.variables[var_name]
            if 'Time' in var.dimensions:
                # we can't support time dependence in the mesh file
                var = time_series_file.variables[var_name]

        else:
            var = time_series_file.variables[var_name]
            
        var_has_time_dim[iVar] ='Time' in var.dimensions

        extra_dim_vals = all_dim_vals[var_name]
        if (extra_dim_vals is None) or (extra_dim_vals.size == 0):
            nHyperSlabs += 1
        else:
            nHyperSlabs += extra_dim_vals.shape[1]

    time_series_file.close()

    any_var_has_time_dim = np.any(var_has_time_dim)

    if any_var_has_time_dim:
        try:
            os.makedirs('%s/time_series'%out_path)
        except OSError:
            pass
    else:
        # there is no point in combining output if no fields have Time dim
        combine_output = False
        nTimes = 1

    # Output time series
    if use_progress_bar:
        widgets = ['Writing time series: ', Percentage(), ' ', Bar(), ' ', ETA()]
        field_bar = ProgressBar(widgets=widgets, maxval=nTimes*nHyperSlabs ).start()
    else:
        print "Writing time series...."

    suffix = 'Cells'
    if any_var_has_time_dim:
        if combine_output or np.all(var_has_time_dim):
            out_prefix = "fieldsOn%s"%suffix
        else:
            out_prefix = "timeDependentFieldsOn%s"%suffix
        # start the pvd file
        pvd_file = utils.write_pvd_header(out_path, out_prefix)
        pvd_file.write('<Collection>\n')

    if not combine_output and not np.all(var_has_time_dim):
        out_prefix = "staticFieldsOn%s"%suffix
        varIndices = np.arange(nVars)[var_has_time_dim == False]
        timeIndependentFile = utils.write_vtp_header(out_path, 
            out_prefix, varIndices[0], varIndices, variable_list,
            all_dim_vals, vertices, connectivity, offsets,
            nPoints, nPolygons, outType, cellData=False, pointData=True)


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
            file_name = '%s/%s.vtp'%(out_path,vtp_file_prefix)
            if append and os.path.exists(file_name):
                continue

            if combine_output:
                varIndices = np.arange(nVars)
            else:
                varIndices = np.arange(nVars)[var_has_time_dim]
            timeDependentFile = utils.write_vtp_header(out_path,
                vtp_file_prefix, varIndices[0], varIndices, variable_list,
                all_dim_vals, vertices, connectivity, offsets,
                nPoints, nPolygons, outType, cellData=False, pointData=True)

            # add time step to pdv file
            pvd_file.write('<DataSet timestep="%d" group="" part="0"\n'%(time_index))
            pvd_file.write('\tfile="%s.vtp"/>\n'%(vtp_file_prefix))

        if time_index == 0 or combine_output:
            varIndices = np.arange(nVars)
        else:
            # only the time-dependent variables
            varIndices = np.arange(nVars)[var_has_time_dim]

        # compute and write zTop and zBot
        (validCellsAtDepth, zIndicesAtDepth, zTop, zBot) = getVertInterpInfo(time_index)

        iHyperSlabProgress = 0
        for iVar in varIndices:
            has_time = var_has_time_dim[iVar]

            var_name = variable_list[iVar]

            if has_time or combine_output:
                vtkFile = timeDependentFile
            else:
                vtkFile = timeIndependentFile
            
            if(var_name == 'zInterfaceTop'):
                vtkFile.appendData(zTop)     
            elif(var_name == 'zInterfaceBot'):
                vtkFile.appendData(zBot)
            else:
            
                (out_var_names, dim_list) = utils.get_hyperslab_name_and_dims(var_name, all_dim_vals)
                for iHyperSlab in range(len(out_var_names)):
                    if dim_list is not None:
                        dim_vals = dim_list[:,iHyperSlab]
                    else:
                        dim_vals = None
    
                    field_var = utils.get_var(var_name, mesh_file, time_series_file)
    
                    field = np.zeros(nCells,dtype=outType)
    
                    try:
                        missing_val = field_var.missing_value
                    except:
                        missing_val = -9999999790214767953607394487959552.000000
    
                    for iBlock in np.arange(0, nBlocks):
                        blockStart = iBlock * blocking
                        blockEnd = min( (iBlock + 1) * blocking, nCells )
                        cellIndices = np.arange(blockStart,blockEnd)
                        field_block = utils.read_field(field_var, dim_vals, 
                                         local_time_indices[time_index], cellIndices,
                                         outType, sign=1, vert_dim_name=vert_dim_name,
                                         nLevels=nLevels)
                                         
                                         
                        if(len(field_block.shape) == 2):
                            valid = validCellsAtDepth[blockStart:blockEnd]
                            zIndices = zIndicesAtDepth[blockStart:blockEnd]
                            fullField = field_block
                            field_block = np.nan*np.ones(blockEnd-blockStart)
                            field_block[valid] = fullField[valid,zIndices[valid]]
    
    
                        field_block[field_block == missing_val] = np.nan
                        field[blockStart:blockEnd] = field_block
                        
                    # map field from cells to points
                    field = field[pointToCellMap]
    
                    vtkFile.appendData(field)
                    
                    del field
                    del field_var

                if use_progress_bar:
                    field_bar.update(time_index*nHyperSlabs + iHyperSlabProgress)
                    iHyperSlabProgress += 1


        if any_var_has_time_dim:
            timeDependentFile.save()
            del timeDependentFile

        if time_index == 0 and not combine_output and not np.all(var_has_time_dim):
            timeIndependentFile.save()
            del timeIndependentFile

    time_series_file.close()
    if use_progress_bar:
        field_bar.finish()

    if any_var_has_time_dim:
        # finish the pdv file
        pvd_file.write('</Collection>\n')
        pvd_file.write('</VTKFile>\n')
#}}}


if __name__ == "__main__":
    if use_progress_bar:
        print " -- Using progress bars --"
    else:
        print " -- Progress bars are not available--"
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file_pattern", dest="filename_pattern", help="MPAS Filename pattern.", metavar="FILE", required=True)
    parser.add_argument("-m", "--mesh_file", dest="mesh_filename", help="MPAS Mesh filename. If not set, it will use the first file in the -f flag as the mesh file.")
    parser.add_argument("--vert_dim", dest="vert_dim_name", help="Vertical dimension", required=True)
    parser.add_argument("--layer_thickness", dest="layer_thickness_name", help="Variable for layer thickness, used to compute zInterface", required=True)
    parser.add_argument("--min_level_cell", dest="min_level_cell_name", help="Index array indicating the minimum valid layer in a cell (default is 0 for all cells)")
    parser.add_argument("--max_level_cell", dest="max_level_cell_name", help="Index array indicating the maximum valid layer in a cell (default is the transect_dim-1 for all cells)")
    parser.add_argument("--z_min", dest="z_min_name", help="Variable specifying the depth of the lower interface (in index space) of the first valid layer on cells, used to compute zInterface")
    parser.add_argument("--z_max", dest="z_max_name", help="Variable specifying the depth of the upper interface (in index space) of the last valid layer on cells, used to compute zInterface")
    parser.add_argument("--field_level_cell", dest="field_level_cell_name", help="Index array on cells indicating the depth at which to extract fields")
    parser.add_argument("--depth", dest="depth", type=float, help="A fixed depth at which to extract fields.")
    parser.add_argument("-b", "--blocking", dest="blocking", default=10000, help="Size of blocks when reading MPAS file", metavar="BLK")
    parser.add_argument("-v", "--variable_list", dest="variable_list", help="List of variables on cells to extract ('all' for all variables on cells)", metavar="VAR", required=True)
    parser.add_argument("-3", "--32bit", dest="output_32bit", help="If set, the vtk files will be written using 32bit floats.", action="store_true")
    parser.add_argument("-c", "--combine", dest="combine_output", help="If set, time-independent fields are written to each file along with time-dependent fields.", action="store_true")
    parser.add_argument("-a", "--append", dest="append", help="If set, only vtp files that do not already exist are written out.", action="store_true")
    parser.add_argument("-d", "--dim_list", dest="dimension_list", nargs="+", help="A list of dimensions and associated values.")
    args = parser.parse_args()

    if not args.output_32bit:
        use_32bit = False
    else:
        use_32bit = True

    out_path = 'vtk_depth_files'
    try:
        os.makedirs(out_path)
    except OSError:
        pass

    (time_indices, time_file_names) = utils.setup_time_indices(args.filename_pattern)

    separate_mesh_file = True
    if not args.mesh_filename:
        args.mesh_filename = time_file_names[0]
        separate_mesh_file = False

    # Setting dimension values:
    time_series_file = NetCDFFile(time_file_names[0], 'r')
    if separate_mesh_file:
        mesh_file = NetCDFFile(args.mesh_filename, 'r')
    else:
        mesh_file = None
    extra_dims = utils.parse_extra_dims( 
        args.dimension_list, time_series_file, mesh_file)
    (all_dim_vals, cellVars, vertexVars, edgeVars) = utils.setup_dimension_values_and_sort_vars(
            time_series_file, mesh_file, args.variable_list, 
            extra_dims, basic_dims=['nCells','Time',args.vert_dim_name],
            include_dims=['nCells'])
    time_series_file.close()
    if(mesh_file is not None):
        mesh_file.close()

    cellVars.append('zInterfaceTop')
    cellVars.append('zInterfaceBot')
    all_dim_vals['zInterfaceTop'] = None
    all_dim_vals['zInterfaceBot'] = None

    utils.summarize_extraction(args.mesh_filename, time_indices, cellVars,
                               vertexVars, edgeVars)

    # Handle cell variables
    if len(cellVars) > 0:
        print " -- Extracting cell fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        # Build cell list
        (vertices, connectivity, offsets, pointToCellMap) = build_point_and_polygon_lists( 
           mesh_file, use_32bit, out_path )

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series( time_indices, time_file_names, mesh_file, args.blocking,
                                 all_dim_vals, cellVars, vertices, connectivity, offsets,
                                 args.vert_dim_name, pointToCellMap, 
                                 args.layer_thickness_name, args.min_level_cell_name,
                                 args.max_level_cell_name,  
                                 args.z_min_name, args.z_max_name, 
                                 args.field_level_cell_name, args.depth,                        
                                 use_32bit, args.combine_output, args.append,
                                 out_path)
        if separate_mesh_file:
            mesh_file.close()

        print ""
        del vertices
        del connectivity
        del offsets


# vim: set expandtab:
