using Test
using ReadVTK, WriteVTK

# Commit in the example file repository for which the test files will be downloaded
# Note: The purpose of using a specific commit hash (instead of `main`) is to be able to tie a given
#       version of ReadVTK to a specific version of the test file repository. This way, also tests
#       for older ReadVTK releases should continue to work.
TEST_EXAMPLES_COMMIT = "92b44ef4666cae5aa5ffe1c5c35b5823c0073c31"

# Local folder to store downloaded example files. If you change this, also adapt `../.gitignore`!
TEST_EXAMPLES_DIR = "examples"


get_test_example_file(filename) = get_example_file(filename, head=TEST_EXAMPLES_COMMIT,
                                                   output_directory=TEST_EXAMPLES_DIR)


# Start with a clean environment: remove example file directory if it exists
isdir(TEST_EXAMPLES_DIR) && rm(TEST_EXAMPLES_DIR, recursive=true)
mkpath(TEST_EXAMPLES_DIR)


@time @testset "ReadVTK" begin
  @testset "basic tests" begin
    @testset "VTKFile" begin
      @test VTKFile(get_test_example_file("celldata_inline_binary_uncompressed.vtu")) isa VTKFile

      mktemp() do path, io
        write(io, "# vtk DataFile Version")
        flush(io)

        @test_throws ErrorException VTKFile(path)
      end
    end

    vtk_file = VTKFile(get_test_example_file("celldata_inline_binary_uncompressed.vtu"))

    @testset "get_cell_data" begin
      @test get_cell_data(vtk_file) isa ReadVTK.VTKData
    end

    cell_data = get_cell_data(vtk_file)

    @testset "VTKData auxiliary functions" begin
      @test firstindex(cell_data) == "cell_ids"
      @test lastindex(cell_data) == "indicator_shock_capturing"
      @test length(cell_data) == 5
      @test size(cell_data) == (5,)
      @test keys(cell_data) == ("cell_ids", "element_ids", "levels",
                                "indicator_amr", "indicator_shock_capturing")
      @test iterate(cell_data) == ("cell_ids" => cell_data["cell_ids"], 2)
      @test_throws KeyError cell_data["does_not_exist"]
      @test eltype(cell_data) == Pair{String, ReadVTK.VTKDataArray}
    end

    @testset "extract VTKDataArray from VTKData" begin
      @test cell_data["cell_ids"] isa ReadVTK.VTKDataArray
    end

    cell_ids = cell_data["cell_ids"]

    @testset "get_data" begin
      @test get_data(cell_ids) isa Base.ReinterpretArray{Int, 1}
    end

    data = get_data(cell_ids)

    @testset "validate data" begin
      @test length(data) == vtk_file.n_cells
      @test first(data) == 4
      @test last(data) == 4113
      @test sum(data) == 6357314
    end

    @testset "get_points" begin
      @test size(get_points(vtk_file)) == (3, 4434)
    end

    points = get_points(vtk_file)

    @testset "validate points" begin
      @test points[:, 1] ≈ [-64.0, -64.0, 0.0]
      @test sum(points) ≈ -6442.5
    end

    @testset "get_cells" begin
      @test get_cells(vtk_file) isa ReadVTK.VTKCells
    end

    cells = get_cells(vtk_file)

    @testset "validate cells" begin
      @test size(cells) == (3085,)
      @test div(Int(sum(cells.types)), 8) == 3085
      @test sum(cells.offsets) == 19040620
      @test cells.connectivity[1000] == 421
    end

    @testset "show" begin
      @test isnothing(show(devnull, vtk_file))
      @test isnothing(show(devnull, cell_data))
      @test isnothing(show(devnull, cell_ids))
    end
  end

  @testset "binary compressed file" begin
    vtk_file = VTKFile(get_test_example_file("celldata_inline_binary_compressed.vtu"))
    cell_data = get_cell_data(vtk_file)
    cell_ids = cell_data["cell_ids"]
    data = get_data(cell_ids)

    @testset "validate data" begin
      @test length(data) == vtk_file.n_cells
      @test first(data) == 4
      @test last(data) == 4113
      @test sum(data) == 6357314
    end
  end

  @testset "appended uncompressed file" begin
    vtk_file = VTKFile(get_test_example_file("celldata_appended_binary_uncompressed.vtu"))
    cell_data = get_cell_data(vtk_file)
    cell_ids = cell_data["cell_ids"]
    data = get_data(cell_ids)

    @testset "validate data" begin
      @test length(data) == vtk_file.n_cells
      @test first(data) == 4
      @test last(data) == 4113
      @test sum(data) == 6357314
    end
  end

  @testset "appended compressed file" begin
    vtk_file = VTKFile(get_test_example_file("celldata_appended_binary_compressed.vtu"))
    cell_data = get_cell_data(vtk_file)
    cell_ids = cell_data["cell_ids"]
    data = get_data(cell_ids)

    @testset "validate data" begin
      @test length(data) == vtk_file.n_cells
      @test first(data) == 4
      @test last(data) == 4113
      @test sum(data) == 6357314
    end
  end

  @testset "point data" begin
    vtk_file = VTKFile(get_test_example_file("pointdata_appended_binary_compressed.vtu"))

    @testset "get_point_data" begin
      @test get_point_data(vtk_file) isa ReadVTK.VTKData
    end

    point_data = get_point_data(vtk_file)
    pressure = point_data["p"]
    data = get_data(pressure)

    @testset "validate data" begin
      @test length(data) == vtk_file.n_points
      @test first(data) ≈ 0.7999332810225936
      @test last(data) ≈ 0.8004962389182811
      @test sum(data) ≈ 192.1204941112099
    end
  end

  # Test for validation of uniform grid ("image data") read feature
  @testset "ImageData" begin
    ## Generate grid file and write vti
    
    # grid geometry parameter
    input_origin  = [ 1.0, 1.0, 2.0 ]
    input_spacing = [ 1.0, 1.0, 0.1 ]
    input_ending  = [ 3.0, 2.0, 2.2 ]

    # compute ranges
    x, y, z = [ (input_origin[i] : input_spacing[i] : input_ending[i])  for i in (1:3) ]
    Nx, Ny, Nz = length(x), length(y), length(z)

    # generate random data
    point_scalar_field = rand(Nx, Ny, Nz)
    point_data_name    = "Point scalar data"

    cell_scalar_field  = rand(Nx-1, Ny-1, Nz-1)
    cell_data_name     = "Name of Test Cell scalar data"

    # write vti file using WriteVTK
    vtk_grid("grid", x, y, z) do vtk
        vtk[ point_data_name, VTKPointData()] = point_scalar_field  # scalar field attached to points
        vtk[ cell_data_name,  VTKCellData() ] = cell_scalar_field   # scalar field attached to cells
    end
    
    ## Read vti file using ReadVTK
    filepath = "grid.vti"
    vtk = VTKFile( filepath )
    
    origin  = get_origin( vtk )
    spacing = get_spacing( vtk )
    cell_data_read = get_data( get_cell_data(vtk)[cell_data_name] )
    point_data_read = get_data( get_point_data(vtk)[point_data_name] )
    
    cell_data_read = reshape( cell_data_read,    ( (Nx-1), (Ny-1), (Nz-1) ) )
    point_data_read = reshape( point_data_read,  (     Nx,     Ny,     Nz ) )

    # test if cell and point data are well read
    @test cell_data_read == cell_scalar_field
    @test point_data_read == point_scalar_field

    # test if the origin and spacing are well read
    @test origin  == input_origin
    @test spacing == input_spacing
  end
end


# Clean up afterwards: delete example file directory
@test_nowarn rm(TEST_EXAMPLES_DIR, recursive=true)
