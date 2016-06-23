#!/usr/bin/env python
"""
Name: paraview_vtk_topo_field_extractor.py
Authors: Xylar Asay-Davis
Date: 04/24/2016

This script extracts geometry for plotting MPAS fields with varying depth 
(presumably topography) as well as fields sampled using a specified index 
field on cells, corresponding to the depth of the topography.

Currently, only fields on cells are supported. In paraview, fields
on cells will have constant values but will be interpolated over
rectangular steps between cells at edges.

Usage is similar to paraview_vtk_field_extractor.py. Additional flags are:
-t, --topo_dim: specifies the vertical dimension of the model (e.g.
  nVertLevels in MPAS-Ocean)

-i, --topo_cell_index: specifies either a constant index value or the name 
  of an index array on cells that corresponds to the location of the 
  topography (e.g. maxLevelCell in MPAS-Ocean). If this flag is excluded, 
  each cell is sampled in the bottom level of topo_dim.

results are stored in topoFieldsOnCells or timeDependentTopoFieldsOnCells
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


def build_point_and_polygon_lists( nc_file, output_32bit ):#{{{
    nCells = len(nc_file.dimensions['nCells'])
    nEdges = len(nc_file.dimensions['nEdges'])

    nEdgesOnCell = nc_file.variables['nEdgesOnCell'][:]
    verticesOnCell = nc_file.variables['verticesOnCell'][:,:]-1
    edgesOnCell = nc_file.variables['edgesOnCell'][:,:]-1
    X_var = nc_file.variables['xVertex']
    Y_var = nc_file.variables['yVertex']
    Z_var = nc_file.variables['zVertex']

    nPoints = np.sum(nEdgesOnCell)
    nPolygons = nCells+nEdges

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
    # a polygon with nEdgesOnCell vertices per cell plus a polygon with 4 vertices per edge
    connectivity = np.zeros(nPoints+nEdges*4, dtype=int)

    # Build cells
    if use_progress_bar:
        widgets = ['Build cell connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        cell_bar = ProgressBar(widgets=widgets, maxval=nCells).start()
    else:
        print "Build cell connectivity..."

    outIndex = 0

    connectivity[0:nPoints] = np.arange(nPoints)

    for iCell in range(nCells):
        nVerts = nEdgesOnCell[iCell]
        iVerts = verticesOnCell[iCell,0:nVerts]
        iEdges = edgesOnCell[iCell,0:nVerts]
        
        X[outIndex:outIndex+nVerts] = X_var[iVerts]
        Y[outIndex:outIndex+nVerts] = Y_var[iVerts]
        Z[outIndex:outIndex+nVerts] = Z_var[iVerts]

        pointToCellMap[outIndex:outIndex+nVerts] = iCell
        for index in range(len(iEdges)):
            iEdge = iEdges[index]
            offset = pointCountOnEdge[iEdge]
            edgeVertToPointMap[iEdges[index], offset] = outIndex + np.mod(index-1, nVerts)
            edgeVertToPointMap[iEdges[index], offset+1] = outIndex + index
            pointCountOnEdge[iEdge] += 2
            

        outIndex += nVerts
        offsets[iCell] = outIndex

        if use_progress_bar:
            cell_bar.update(iCell)

    if use_progress_bar:
        cell_bar.finish()

    # Build cells
    if use_progress_bar:
        widgets = ['Build edge connectivity: ', Percentage(), ' ', Bar(), ' ', ETA()]
        cell_bar = ProgressBar(widgets=widgets, maxval=nEdges).start()
    else:
        print "Build cell connectivity..."

    for iEdge in range(nEdges):
        connectivity[outIndex] = edgeVertToPointMap[iEdge,0]
        connectivity[outIndex+1] = edgeVertToPointMap[iEdge,1]
        if(pointCountOnEdge[iEdge] > 2):
            connectivity[outIndex+2] = edgeVertToPointMap[iEdge,2]
            connectivity[outIndex+3] = edgeVertToPointMap[iEdge,3]
        else:
            connectivity[outIndex+2] = edgeVertToPointMap[iEdge,1]
            connectivity[outIndex+3] = edgeVertToPointMap[iEdge,0]
        
        outIndex += 4
        offsets[nCells + iEdge] = outIndex

        if use_progress_bar:
            cell_bar.update(iEdge)

    if use_progress_bar:
        cell_bar.finish()


    return ((X,Y,Z), connectivity, offsets, pointToCellMap)

#}}}


def build_field_time_series( local_time_indices, file_names, mesh_file, all_dim_vals,
                             variable_list, vertices, connectivity, offsets,
                             topo_dim, topo_cell_indices, pointToCellMap, 
                             output_32bit, combine_output, append ):#{{{

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
        
    if (mesh_file is not None) and (topo_dim in mesh_file.dimensions):
        nTopoLevels = len(mesh_file.dimensions[topo_dim])
    else:
        nTopoLevels = len(time_series_file.dimensions[topo_dim])

    nPolygons = len(offsets)
    nPoints = len(vertices[0])
    nTimes = len(local_time_indices)
    nVars = len(variable_list)
    
    var_has_time_dim = np.zeros(nVars,bool)
    nHyperSlabs = 0
    for iVar in range(nVars):
        var_name = variable_list[iVar]
        if (mesh_file is not None) and (var_name in mesh_file.variables):
            # we can't support time dependence in the mesh file
            assert('Time' not in mesh_file.variables[var_name].dimensions)
            var_has_time_dim[iVar] = False
        else:
            var_has_time_dim[iVar] = 'Time' in time_series_file.variables[var_name].dimensions

        extra_dim_vals = all_dim_vals[var_name]
        if (extra_dim_vals is None) or (extra_dim_vals.size == 0):
            nHyperSlabs += 1
        else:
            nHyperSlabs += extra_dim_vals.shape[1]

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
        field_bar = ProgressBar(widgets=widgets, maxval=nTimes*nHyperSlabs ).start()
    else:
        print "Writing time series...."

    suffix = 'Cells'
    if any_var_has_time_dim:
        if combine_output or np.all(var_has_time_dim):
            out_prefix = "topoFieldsOn%s"%suffix
        else:
            out_prefix = "timeDependentTopoFieldsOn%s"%suffix
        # start the pvd file
        pvd_file = utils.write_pvd_header(out_prefix)
        pvd_file.write('<Collection>\n')

    if not combine_output and not np.all(var_has_time_dim):
        out_prefix = "staticTopoFieldsOn%s"%suffix
        varIndices = np.arange(nVars)[var_has_time_dim == False]
        timeIndependentFile = utils.write_vtp_header(
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
            file_name = 'vtk_files/%s.vtp'%(vtp_file_prefix)
            if append and os.path.exists(file_name):
                continue

            if combine_output:
                varIndices = np.arange(nVars)
            else:
                varIndices = np.arange(nVars)[var_has_time_dim]
            timeDependentFile = utils.write_vtp_header(
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


        iHyperSlabProgress = 0
        for iVar in varIndices:
            use_time_dependent = var_has_time_dim[iVar] or combine_output

            var_name = variable_list[iVar]
            if (mesh_file is not None) and (var_name in mesh_file.variables):
                nc_file = mesh_file
            else:
                nc_file = time_series_file
            field_var = nc_file.variables[var_name]

            field_ndims = len(field_var.dimensions)
            

            (out_var_names, dim_list) = utils.get_hyperslab_name_and_dims(
                var_name, all_dim_vals)                    
                        
            if use_time_dependent:
                vtkFile = timeDependentFile
            else:
                vtkFile = timeIndependentFile
            
            for iHyperSlab in range(len(out_var_names)):
                
                try:
                    missing_val = field_var.missing_value
                except:
                    missing_val = -9999999790214767953607394487959552.000000

                # Build data
                dim_vals = []
                extraDimIndex = 0
                for iDim in range(field_ndims):
                    dim = field_var.dimensions[iDim]
                    if dim == 'Time':
                        dim_vals.append(local_time_indices[time_index])
                    elif dim == 'nCells':
                        dim_vals.append(np.arange(nCells))
                    elif dim == topo_dim:
                        dim_vals.append(np.arange(nTopoLevels))
                    else:
                        dim_vals.append(dim_list[extraDimIndex,iHyperSlab])
                        extraDimIndex += 1
                        
                if field_ndims == 1:
                    field = field_var[dim_vals[0]]
                elif field_ndims == 2:
                    field = field_var[dim_vals[0], dim_vals[1]]
                elif field_ndims == 3:
                    field = field_var[dim_vals[0], dim_vals[1], dim_vals[2]]
                elif field_ndims == 4:
                    field = field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3]]
                elif field_ndims == 5:
                    field = field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3], dim_vals[4]]

                if topo_dim in field_var.dimensions:
                    field = field[np.arange(nCells),topo_cell_indices]
                
                # convert to the same type as field before masking with NaNs
                field = np.array(field, dtype=outType)
                field[field == missing_val] = np.nan
                
                # map field from cells to points
                field = field[pointToCellMap]

                vtkFile.appendData(field)

                if use_progress_bar:
                    field_bar.update(time_index*nHyperSlabs + iHyperSlabProgress)
                    iHyperSlabProgress += 1

                del field
                del field_ndims
                del field_var

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
    parser.add_argument("-t", "--topo_dim", dest="topo_dim", help="Vertical dimension", required=True)
    parser.add_argument("-i", "--topo_cell_index", dest="topo_cell_index", help="Index array indicating the depth of the topograph (default is the topo_dim-1 for all cells)")
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
    (topo_cell_indices, extra_dims) = utils.parse_extra_dims( 
        args.dimension_list, time_series_file, mesh_file,
        topo_dim=args.topo_dim, topo_cell_indices_name=args.topo_cell_index)
    (all_dim_vals, cellVars, vertexVars, edgeVars) = utils.setup_dimension_values_and_sort_vars(
            time_series_file, mesh_file, args.variable_list, 
            extra_dims, basic_dims=['nCells','Time',args.topo_dim],
            include_dims=['nCells'])
    time_series_file.close()
    if(mesh_file is not None):
        mesh_file.close()

    utils.summarize_extraction(args.mesh_filename, time_indices, cellVars,
                               vertexVars, edgeVars)

    # Handle cell variables
    if len(cellVars) > 0:
        print " -- Extracting cell fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        # Build cell list
        (vertices, connectivity, offsets, pointToCellMap) = build_point_and_polygon_lists( 
           mesh_file, use_32bit )

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series( time_indices, time_file_names, mesh_file,
                                 all_dim_vals, cellVars, vertices, connectivity, offsets,
                                 args.topo_dim, topo_cell_indices, pointToCellMap, 
                                 use_32bit, args.combine_output, args.append )
        if separate_mesh_file:
            mesh_file.close()

        print ""
        del vertices
        del connectivity
        del offsets


# vim: set expandtab:
