#!/usr/bin/env python
"""
Name: paraview_vtk_field_extractor.py
Authors: Xylar Asay-Davis
Date: 06/20/2016

"""
import os
import numpy as np

from netCDF4 import Dataset as NetCDFFile
import argparse

try:
    from progressbar import ProgressBar, Percentage, Bar, ETA
    use_progress_bar = True
except:
    use_progress_bar = False

import utils


def build_field_time_series( local_time_indices, file_names, mesh_file, blocking, all_dim_vals,
                             blockDimName, variable_list, vertices, connectivity, offsets,
                             vert_dim_name, layer_thickness_name,
                             min_level_name, max_level_name,
                             z_min_name, z_max_name, depth,
                             valid_mask, output_32bit, combine_output, append ):#{{{
                             
                             
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
            
        cellIndices = np.arange(blockDim)
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
            zMax = np.zeros(blockDim)

        zInterface = utils.compute_zInterface(minLevelCell, maxLevelCell, layerThickness, 
                                              zMin, zMax, dtype=outType)
        valid = np.zeros(blockDim,bool)
        zIndices = -1*np.ones(blockDim,int)
        for zIndex in range(nLevels):
            if(layerSign == 1):
                mask = np.logical_and(depth >= zInterface[:,zIndex], 
                                      depth < zInterface[:,zIndex+1])
            else:
                mask = np.logical_and(depth <= zInterface[:,zIndex], 
                                      depth > zInterface[:,zIndex+1])
            valid[mask] = True
            zIndices[mask] = zIndex
            
        return (valid, zIndices)

        

    if len(variable_list) == 0:
        return

    if output_32bit:
        outType = 'float32'
    else:
        outType = 'float64'


    # Get dimension info to allocate the size of Colors
    time_series_file = NetCDFFile(file_names[0], 'r')

    if mesh_file is not None:
        # blockDim may not exist in time series file
        blockDim = len(mesh_file.dimensions[blockDimName])
    else:
        blockDim = len(time_series_file.dimensions[blockDimName])

    # Pre-compute the number of blocks
    nBlocks = 1 + blockDim / blocking


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
        
        
    (validCellsAtDepth, zIndicesAtDepth) = getVertInterpInfo(time_index=0)
        
    var_has_time_dim = np.zeros(nVars,bool)
    nHyperSlabs = 0
    for iVar in range(nVars):
        var_name = variable_list[iVar]
        if var_name in time_series_file.variables:
            var_has_time_dim[iVar] = 'Time' in time_series_file.variables[var_name].dimensions
        else:
            # we can't support time dependence in the mesh file
            assert('Time' not in mesh_file.variables[var_name].dimensions)
            var_has_time_dim[iVar] = False

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

    suffix = blockDimName[1:]
    if any_var_has_time_dim:
        if combine_output or np.all(var_has_time_dim):
            out_prefix = "fieldsOn%s"%suffix
        else:
            out_prefix = "timeDependentFieldsOn%s"%suffix
        # start the pvd file
        pvd_file = utils.write_pvd_header(out_prefix)
        pvd_file.write('<Collection>\n')

    if not combine_output and not np.all(var_has_time_dim):
        out_prefix = "staticFieldsOn%s"%suffix
        varIndices = np.arange(nVars)[var_has_time_dim == False]
        timeIndependentFile = utils.write_vtp_header(out_prefix, varIndices[0], varIndices,
                                                     variable_list, all_dim_vals,
                                                     vertices, connectivity, offsets, 
                                                     nPoints, nPolygons, outType, 
                                                     cellData=True, pointData=False)


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
            timeDependentFile = utils.write_vtp_header(vtp_file_prefix, varIndices[0], varIndices,
                                                       variable_list, all_dim_vals,
                                                       vertices, connectivity, offsets, 
                                                       nPoints, nPolygons, outType, 
                                                       cellData=True, pointData=False)

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
            has_time = var_has_time_dim[iVar]

            var_name = variable_list[iVar]
            (out_var_names, dim_list) = utils.get_hyperslab_name_and_dims(var_name, all_dim_vals)
            if has_time or combine_output:
                vtkFile = timeDependentFile
            else:
                vtkFile = timeIndependentFile
            for iHyperSlab in range(len(out_var_names)):
                if dim_list is not None:
                    dim_vals = dim_list[:,iHyperSlab]
                else:
                    dim_vals = None

                field_var = utils.get_var(var_name, mesh_file, time_series_file)

                field = np.zeros(blockDim,dtype=outType)

                try:
                    missing_val = field_var.missing_value
                except:
                    missing_val = -9999999790214767953607394487959552.000000

                for iBlock in np.arange(0, nBlocks):
                    blockStart = iBlock * blocking
                    blockEnd = min( (iBlock + 1) * blocking, blockDim )
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

                field = field[valid_mask]

                vtkFile.appendData(field)

                if use_progress_bar:
                    field_bar.update(time_index*nHyperSlabs + iHyperSlabProgress)
                    iHyperSlabProgress += 1

                del field
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
    parser.add_argument("--vert_dim", dest="vert_dim_name", help="Vertical dimension", required=True)
    parser.add_argument("--layer_thickness", dest="layer_thickness_name", help="Variable for layer thickness, used to compute zInterface", required=True)
    parser.add_argument("--min_level_cell", dest="min_level_cell_name", help="Index array indicating the minimum valid layer in a cell (default is 0 for all cells)")
    parser.add_argument("--max_level_cell", dest="max_level_cell_name", help="Index array indicating the maximum valid layer in a cell (default is the transect_dim-1 for all cells)")
    parser.add_argument("--z_min", dest="z_min_name", help="Variable specifying the depth of the lower interface (in index space) of the first valid layer on cells, used to compute zInterface")
    parser.add_argument("--z_max", dest="z_max_name", help="Variable specifying the depth of the upper interface (in index space) of the last valid layer on cells, used to compute zInterface")
    parser.add_argument("--depth", dest="depth", type=float, default=0.0, help="Depth at which to extract fields.")
    parser.add_argument("-b", "--blocking", dest="blocking", help="Size of blocks when reading MPAS file", metavar="BLK")
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

    if not args.blocking:
        args.blocking = int(10000)
    else:
        args.blocking = int(args.blocking)

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
    extra_dims = utils.parse_extra_dims(args.dimension_list, time_series_file, mesh_file)
    # only allow nCells dims; exclude vert_dim_name from extra_dims
    (all_dim_vals, cellVars, vertexVars, edgeVars) = utils.setup_dimension_values_and_sort_vars(
            time_series_file, mesh_file, args.variable_list, extra_dims,
            basic_dims=['nCells', 'nEdges', 'nVertices', 'Time', args.vert_dim_name],
            include_dims=['nCells'])
    time_series_file.close()
    if(mesh_file is not None):
        mesh_file.close()

    utils.summarize_extraction(args.mesh_filename, time_indices, cellVars, vertexVars, edgeVars)

    # Handle cell variables
    if len(cellVars) > 0:
        print " -- Extracting cell fields --"

        mesh_file = NetCDFFile(args.mesh_filename, 'r')

        # Build vertex list
        vertices = utils.build_location_list_xyz( mesh_file, 'xVertex', 'yVertex', 'zVertex', use_32bit )

        # Build cell list
        (connectivity, offsets, valid_mask) = utils.build_cell_lists( mesh_file, args.blocking )

        if not separate_mesh_file:
            mesh_file.close()
            mesh_file = None

        build_field_time_series( time_indices, time_file_names, mesh_file, args.blocking,
                                 all_dim_vals, 'nCells', cellVars, vertices, connectivity, offsets,
                                 args.vert_dim_name, args.layer_thickness_name,
                                 args.min_level_cell_name, args.max_level_cell_name,
                                 args.z_min_name, args.z_max_name, args.depth,
                                 valid_mask, use_32bit, args.combine_output, args.append )
        if separate_mesh_file:
            mesh_file.close()


# vim: set expandtab:
