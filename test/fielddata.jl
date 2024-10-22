# Create some cells
points = rand(3, 5)
cells = [
     MeshCell(VTKCellTypes.VTK_TRIANGLE, [1, 4, 2]),
     MeshCell(VTKCellTypes.VTK_QUAD, [2, 4, 3, 5]),
   ]

scalar = 5.
vector = [1., 2., 3.]

path = joinpath(TEST_EXAMPLES_DIR, "fields")

vtk_grid(path, points, cells, append=false) do vtk
    vtk["scalar"] = scalar
    vtk["vector"] = vector
end

vtk_read = VTKFile(path * ".vtu")
fielddata = get_field_data(vtk_read)

@test only(get_data(fielddata["scalar"])) == scalar

@test get_data(fielddata["vector"]) == vector