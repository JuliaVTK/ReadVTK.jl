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

function create_directory(TEST_EXAMPLES_DIR)
  isdir(TEST_EXAMPLES_DIR) && rm(TEST_EXAMPLES_DIR, recursive=true)
  mkpath(TEST_EXAMPLES_DIR)
  return nothing
end                                                   

clean_directory(TEST_EXAMPLES_DIR) = @test_nowarn rm(TEST_EXAMPLES_DIR, recursive=true)


@time @testset "ReadVTK" begin
  @testset "basic tests" begin
    # Start with a clean environment: remove example file directory if it exists
    create_directory(TEST_EXAMPLES_DIR)
    
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

    # Clean up afterwards: delete example file directory
    clean_directory(TEST_EXAMPLES_DIR)
  end

  # Test for validation of uniform grid ("image data") read feature
  @testset "ImageData" begin
    # Start with a clean environment: remove example file directory if it exists
    create_directory(TEST_EXAMPLES_DIR)

    ## Generate grid file and write vti
    
    # grid geometry parameters
    input_origin   = [1.0, 1.0, 2.0]
    input_ending   = [3.0, 2.0, 2.2]
    input_numNodes = [  4,   2,   2]  

    # compute ranges
    input_spacing = (input_ending .- input_origin) ./ (input_numNodes.-1)
    x, y, z = [(input_origin[i]:input_spacing[i]:input_ending[i]) for i in (1:3)]
    Nx, Ny, Nz = length(x), length(y), length(z)

    # generate random data
    point_scalar_field = rand(Nx, Ny, Nz)
    point_data_name    = "Point scalar data"

    cell_scalar_field  = rand(Nx-1, Ny-1, Nz-1)
    cell_data_name     = "Cell scalar data"

    # write vti file using WriteVTK
    path = joinpath(TEST_EXAMPLES_DIR, "grid")
    vtk_grid(path, x, y, z) do vtk
        vtk[point_data_name, VTKPointData()] = point_scalar_field # scalar field attached to points
        vtk[cell_data_name, VTKCellData()] = cell_scalar_field    # scalar field attached to cells
    end
    
    # Read vti file using ReadVTK
    filepath = joinpath(TEST_EXAMPLES_DIR, "grid.vti")
    @testset "VTKFile" begin
      @test VTKFile(filepath) isa VTKFile
    end
    vtk = VTKFile(filepath)

    # check getter functions
    @testset "get_origin" begin
      @test get_origin(vtk) == input_origin
    end

    @testset "get_spacing" begin
      @test get_spacing(vtk) == input_spacing
    end

    @testset "get scalar cell data" begin
      cell_data_raw = get_data(get_cell_data(vtk)[cell_data_name])
      cell_data_reshaped = reshape(cell_data_raw, ((Nx-1), (Ny-1), (Nz-1)))
      cell_data_reshaped1 = get_data_reshaped(get_cell_data(vtk)[cell_data_name], cell_data=true)
      
      @test cell_data_reshaped == cell_scalar_field
      @test cell_data_reshaped1 == cell_scalar_field
    end

    @testset "get scalar point data" begin
      point_data_raw = get_data(get_point_data(vtk)[point_data_name])
      point_data_reshaped = reshape(point_data_raw, (Nx, Ny, Nz))
      point_data_reshaped1 =  get_data_reshaped(get_point_data(vtk)[point_data_name])
      
      @test point_data_reshaped == point_scalar_field
      @test point_data_reshaped1 == point_scalar_field
    end

    # generate random 2D data
    point_scalar_field = rand(Nx, Ny)
    cell_scalar_field  = rand(Nx-1, Ny-1)
    
    # write 2D vti file using WriteVTK
    path = joinpath(TEST_EXAMPLES_DIR, "grid_2D")
    vtk_grid(path, x, y) do vtk
      vtk[point_data_name, VTKPointData()] = point_scalar_field # scalar field attached to points
      vtk[cell_data_name, VTKCellData()] = cell_scalar_field    # scalar field attached to cells
    end

    filepath = joinpath(TEST_EXAMPLES_DIR, "grid_2D.vti")
    @testset "VTKFile2D" begin
      @test VTKFile(filepath) isa VTKFile
    end
    vtk = VTKFile(filepath)

    @testset "get_2D_origin" begin
      @test get_origin(vtk) == input_origin[1:2]
    end

    @testset "get_2D_spacing" begin
      @test get_spacing(vtk) == input_spacing[1:2]
    end

    @testset "get 2D scalar cell data" begin
      cell_data_raw = get_data(get_cell_data(vtk)[cell_data_name])
      cell_data_reshaped = reshape(cell_data_raw, ((Nx-1), (Ny-1)))
      @test cell_data_reshaped == cell_scalar_field
    end

    @testset "get 2D scalar point data" begin
      point_data_raw = get_data(get_point_data(vtk)[point_data_name])
      point_data_reshaped = reshape(point_data_raw, (Nx, Ny))
      @test point_data_reshaped == point_scalar_field
    end

    # Clean up afterwards: delete example file directory
    clean_directory(TEST_EXAMPLES_DIR)
  end

  # Test set for PolyData
  @testset "PolyData" begin
    # Start with a clean environment: remove example file directory if it exists
    create_directory(TEST_EXAMPLES_DIR)

    # function to convert VTKPrimitives to a sequence of Int arrays
    function primitives_to_arrays(primitives::VTKPrimitives)::Vector
      out = []
      first = 1
      for last in primitives.offsets
        push!(out, primitives.connectivity[first:last] .+ 1)
        first = last + 1
      end
      return out
    end

    @testset "single type" begin
      
      # define some points
      n = 20
      points = zeros(3, n)
      points[1, :] .= [cos(4 * pi * i / n) for i in 0:(n - 1)]
      points[2, :] .= [sin(4 * pi * i / n) for i in 0:(n - 1)]
      points[3, :] .= [0.2 * i for i in 0:n-1]

      # define polygons
      polys = [i:(i + 3) for i = 1:(n - 3)]
      cells = [MeshCell(PolyData.Polys(), p) for p in polys]

      # define values
      point_values = [4 * pi * i / n for i in 0:(n - 1)]
      cell_values = [0.3 * i for i in 1:(n - 3)]

      # save using WriteVTK
      path = joinpath(TEST_EXAMPLES_DIR, "spiral")
      vtk_grid(path, points, cells) do vtk
        vtk["theta", VTKPointData()] = point_values # scalar field attached to points
        vtk["h", VTKCellData()] = cell_values       # scalar field attached to cells
      end

      # read data from the vtp file
      vtk = VTKFile(path * ".vtp")

      # test correctness
      @test points == get_points(vtk)
      @test point_values == get_data(get_point_data(vtk)["theta"])
      @test cell_values == get_data(get_cell_data(vtk)["h"])
      @test polys == primitives_to_arrays(get_primitives(vtk, "Polys"))
      @test_throws Exception get_primitives(vtk, "Foo")
      @test_throws Exception get_primitives(vtk, "Verts")
    end

    @testset "mixed types" begin
      # define points of a regular tetrahedron
      isqrt2 = 1 / sqrt(2)
      points = permutedims(Float32[
          1  0 -isqrt2 #1
         -1  0 -isqrt2 #2
          0 -1  isqrt2 #3
          0  1  isqrt2 #4
      ])

      # define verts
      verts = [[1], [2], [3], [4]]
      cells_verts = [MeshCell(PolyData.Verts(), v) for v in verts]

      # define lines
      lines = [[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [4, 1]]
      cells_lines = [MeshCell(PolyData.Lines(), l) for l in lines]
      
      # define polys
      polys = [[1, 2, 3], [2, 3, 4], [3, 4, 1], [4, 1, 2]]
      cells_polys = [MeshCell(PolyData.Polys(), p) for p in polys]

      # define cell values
      cell_values = [i for i in 1:14]

      # save using WriteVTK
      path = joinpath(TEST_EXAMPLES_DIR, "tetrahedron")
      vtk_grid(path, points, cells_verts, cells_lines, cells_polys) do vtk
        vtk["id", VTKCellData()] = cell_values
      end

      # read data from the vtp file
      vtk = VTKFile(path*".vtp")

      # test correctness
      @test cell_values == get_data(get_cell_data(vtk)["id"])
      @test verts == primitives_to_arrays(get_primitives(vtk, "Verts"))
      @test lines == primitives_to_arrays(get_primitives(vtk, "Lines"))
      @test polys == primitives_to_arrays(get_primitives(vtk, "Polys"))
    
    end
    
    # Clean up afterwards: delete example file directory
    clean_directory(TEST_EXAMPLES_DIR)
  end 
  
  @testset "RectilinearGrid" begin
    # Start with a clean environment: remove example file directory if it exists
    create_directory(TEST_EXAMPLES_DIR)

    include("rectilinear.jl")

    # Clean up afterwards: delete example file directory
    clean_directory(TEST_EXAMPLES_DIR)
  end

  @testset "Parallel VTK (PVTK) files" begin
    # Start with a clean environment: remove example file directory if it exists
    create_directory(TEST_EXAMPLES_DIR)

    include("pvtk_files.jl")

    # Clean up afterwards: delete example file directory
    clean_directory(TEST_EXAMPLES_DIR)
  end
  
  @testset "StructuredGrid" begin
    # Start with a clean environment: remove example file directory if it exists
    create_directory(TEST_EXAMPLES_DIR)

    include("structuredgrid.jl")

    # Clean up afterwards: delete example file directory
    clean_directory(TEST_EXAMPLES_DIR)
  end

end