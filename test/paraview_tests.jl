
vtk_case1 = VTKFile(get_test_example_file("celldata_inline_binary_uncompressed.vtu"))
vtk_case2 = VTKFile(get_test_example_file("celldata_inline_binary_uncompressed_ParaView.vtu"))

points1 = get_points(vtk_case1)
points2 = get_points(vtk_case2)

cells1 = get_cells(vtk_case1)
cells2 = get_cells(vtk_case2)

cell_data1 = get_cell_data(vtk_case1)
cell_data2 = get_cell_data(vtk_case2)

vtk_keys1 = keys(cell_data1)
vtk_keys2 = keys(cell_data2)

@testset "read points and cells" begin

  @test points1 == points2
  @test cells1.connectivity == cells2.connectivity
  @test cells1.offsets == cells2.offsets
  @test cells1.types == cells2.types

end

@testset "read cell data" begin

  @test vtk_keys1 == vtk_keys2

  for (key1, key2) in zip(vtk_keys1, vtk_keys2)
    @test get_data(cell_data1[key1]) == get_data(cell_data2[key2])
  end

end

write_vtk_out1 = joinpath(TEST_EXAMPLES_DIR, "write_vtk_out1")
write_vtk_out2 = joinpath(TEST_EXAMPLES_DIR, "write_vtk_out2")

vtk_grid(write_vtk_out1, points1, to_meshcells(cells1)) do vtk

  for key1 in vtk_keys1
    vtk[key1] = get_data(cell_data1[key1])
  end

end

vtk_grid(write_vtk_out2, points2, to_meshcells(cells2)) do vtk

  for key2 in vtk_keys2
    vtk[key2] = get_data(cell_data2[key2])
  end

end

write_vtk_out1_UInt8 = read(write_vtk_out1 * ".vtu")
write_vtk_out2_UInt8 = read(write_vtk_out2 * ".vtu")

@testset "write consistency" begin

  @test write_vtk_out1_UInt8 == write_vtk_out2_UInt8

end
