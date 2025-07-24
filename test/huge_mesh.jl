points = rand(3, 1000_000)
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [i, i + 1, i + 2]) for i in 1:999998]

path = joinpath(TEST_EXAMPLES_DIR, "huge_mesh")

vtk_grid(path, points, cells, append = false) do vtk
end

@test_nowarn vtk_read = VTKFile(path * ".vtu")
