
# ReadVTK.jl gridExample
using WriteVTK, ReadVTK

## Generate grid file and write vti
x, y, z = 1:3, 1:2, 2:0.1:2.2
Nx, Ny, Nz = length(x), length(y), length(z)

pointScalarField = rand(Nx, Ny, Nz)
cellScalarField  = rand(Nx-1, Ny-1, Nz-1)

print( "cell scalar field", cellScalarField,"\n")

cellDataName  = "Name of Test Cell scalar data"
pointDataName = "Point scalar data"

print("  writing vtk file...")
vtk_grid("grid", x, y, z) do vtk
    vtk[ pointDataName, VTKPointData()] = pointScalarField  # scalar field attached to points
    vtk[ cellDataName, VTKCellData()] = cellScalarField     # scalar field attached to cells
end
print("done.\n")

## Read vti file
print("  reading vtk file...")
filepath = "grid.vti"
vtk = VTKFile( filepath )
data = get_data( get_cell_data(vtk)[cellDataName] )
print("done.\n")

reshapedData = reshape( cellScalarField, ( (Nx-1), (Ny-1), (Nz-1) ) )

difference = reshapedData .- cellScalarField

