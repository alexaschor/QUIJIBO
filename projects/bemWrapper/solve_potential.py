#!/usr/bin/env python3
# coding: utf-8

import sys
import bempp.api
import numpy as np
import struct

def readDimensionsFromFile(filename):
    with open(filename, 'rb') as f:

        # FIXME Assumes sizeof(int)==4, sizeof(double)==8, w/ little-endian byte order
        # As currently written, the sizes depend on the compiler and machine that ran the
        # previous steps. This could create issues, particularly because of the Python/C++
        # interface.

        xRes, = struct.unpack('<i', f.read(4))
        yRes, = struct.unpack('<i', f.read(4))
        zRes, = struct.unpack('<i', f.read(4))

        xCenter, = struct.unpack('<d', f.read(8))
        yCenter, = struct.unpack('<d', f.read(8))
        zCenter, = struct.unpack('<d', f.read(8))

        xLength, = struct.unpack('<d', f.read(8))
        yLength, = struct.unpack('<d', f.read(8))
        zLength, = struct.unpack('<d', f.read(8))

        print(" Field res:     ", xRes, yRes, zRes)
        print(" Field center:  ", xCenter, yCenter, zCenter)
        print(" Field lengths: ", xLength, yLength, zLength)

        return ((xRes, yRes, zRes),
                (xCenter, yCenter, zCenter),
                (xLength, yLength, zLength))

# Dirichlet trace data function
@bempp.api.real_callable
def _dirichlet_data_trace(x, n, domain_index, result):
    result[0] = 1


def getPotentialField(mesh, dimensions):
    # NOTE: The original C++ implementation (with BEM++ 2.0) set some additional parameters
    # which seem to no longer exist (at least not under the same names) in Bempp-CL. The
    # program seems to work just fine without them, but the parameters omitted here were:

    #  accuracyOptions.doubleRegular.setRelativeQuadratureOrder(8);
    #  accuracyOptions.singleRegular.setRelativeQuadratureOrder(8);
    #
    #  acaOptions.eps = 1e-6;
    #  acaOptions.maximumBlockSize = (int)(grid->leafView()->entityCount(0) / 8);


    # Initialize the spaces
    pwiseLinears = bempp.api.function_space(mesh, "P", 1)
    pwiseConstants = bempp.api.function_space(mesh, "DP", 0)

    # Construct elementary operators
    slpOp = bempp.api.operators.boundary.laplace.single_layer(pwiseConstants, pwiseLinears, pwiseConstants)
    dlpOp = bempp.api.operators.boundary.laplace.double_layer(pwiseLinears, pwiseLinears, pwiseConstants)
    idOp  = bempp.api.operators.boundary.sparse.identity(pwiseLinears, pwiseLinears, pwiseConstants)

    # Construct the grid function representing the (input) Dirichlet data
    dirichlet_trace = bempp.api.GridFunction(pwiseLinears, fun=_dirichlet_data_trace)

    # Construct the right-hand-side grid function
    print(" Constructing RHS grid function...")
    rhs = (-0.5 * idOp + dlpOp) * dirichlet_trace

    # Solve the equation
    print(" Solving equation using GMRES...")
    solFun, info = bempp.api.linalg.gmres(slpOp, rhs, tol=1e-08)

    # Export as VTK
    print(" Exporting solution as grid function...")
    bempp.api.export("solution.vtk", grid_function=solFun, data_type='element')

    # Now that we have a solution, we need to generate a potential field to
    # export to the next step in the pipeline
    print(" Generating 3D potential field...")

    # Unpack the dimensions from the input
    ((xRes,    yRes,    zRes   ),
     (xCenter, yCenter, zCenter),
     (xLength, yLength, zLength)) = dimensions

    xHalfLength = 0.5 * xLength
    yHalfLength = 0.5 * yLength
    zHalfLength = 0.5 * zLength

    dx = xLength / xRes
    dy = yLength / yRes
    dz = zLength / zRes

    # Create list of points
    evaluationPoints3D = np.zeros((3, xRes * yRes * zRes))
    for z in range(0, zRes):
        for y in range(0, yRes):
            for x in range(0, xRes):
                index = x + y * xRes + z * xRes * yRes;

                xReal = xCenter - xHalfLength + x * dx + 0.5 * dx
                yReal = yCenter - yHalfLength + y * dy + 0.5 * dy
                zReal = zCenter - zHalfLength + z * dz + 0.5 * dz

                evaluationPoints3D[0, index] = xReal
                evaluationPoints3D[1, index] = yReal
                evaluationPoints3D[2, index] = zReal

    # Now we need to evaluate each of those points
    print(" Evaluating", xRes * yRes * zRes, "points in potential field...")

    # Create potential operators
    slPotOp = bempp.api.operators.potential.laplace.single_layer(pwiseConstants, evaluationPoints3D)
    dlPotOp = bempp.api.operators.potential.laplace.double_layer(pwiseLinears, evaluationPoints3D)

    # I believe this is Green's representation theorem?
    field3D = dlPotOp * dirichlet_trace - slPotOp * solFun

    return field3D


def writeField3D(filename, field3D, dimensions):
    # Unpack the dimensions
    ((xRes,    yRes,    zRes   ),
     (xCenter, yCenter, zCenter),
     (xLength, yLength, zLength)) = dimensions

    print(f" Writing field to {filename}...")

    with open(filename, 'wb') as f:
        f.write(struct.pack('<i', xRes))
        f.write(struct.pack('<i', yRes))
        f.write(struct.pack('<i', zRes))

        f.write(struct.pack('<d', xCenter))
        f.write(struct.pack('<d', yCenter))
        f.write(struct.pack('<d', zCenter))

        f.write(struct.pack('<d', xLength))
        f.write(struct.pack('<d', yLength))
        f.write(struct.pack('<d', zLength))

        totalCells = xRes * yRes * zRes

        for i in range(totalCells):
            f.write(struct.pack('<d', 1.0 - field3D[0, i]))


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print(" USAGE:", sys.argv[0], " <input Gmsh filename> <input distance field (just for the dimensions)> <output field3D filename>");
        exit(1)

    meshFilename = sys.argv[1]
    field3DFilename = sys.argv[2]
    outputFilename = sys.argv[3]

    print(" Using the Gmsh file:", meshFilename)
    print(" Using the dimensions in the FIELD_3D file:", field3DFilename)
    print(" Outputting potential field to:", outputFilename)

    # Import the mesh file
    print(f" Reading mesh from file...")
    mesh = bempp.api.import_grid(meshFilename)

    # Import the field3D just to read dimensions
    dimensions = readDimensionsFromFile(field3DFilename)

    # Solve for the potential, returns a grid function
    potentialField = getPotentialField(mesh, dimensions)

    # Write the potential field to a file
    writeField3D(outputFilename, potentialField, dimensions)

    print(" Done.")
