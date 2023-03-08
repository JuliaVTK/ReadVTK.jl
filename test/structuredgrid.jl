using StaticArrays: SVector, SMatrix

const FloatType = Float32
const vtk_filename_struct_noext = joinpath(TEST_EXAMPLES_DIR, "structuredgrid")

outfiles = String[]

for compress in [true, false]
  for dim in 2:3
    # Define grid.
    if dim == 2
      Ni, Nj, Nk = 20, 30, 1

      x = zeros(FloatType, Ni,Nj)
      y = zeros(FloatType, Ni,Nj)
  
      for i = 1:Ni, j = 1:Nj
        x[i,j] = i*i/Ni/Ni
        y[i,j] = sqrt(j/Nj)
      end
    elseif dim == 3
      Ni, Nj, Nk = 20, 30, 40
      x = zeros(FloatType, Ni,Nj,Nk)
      y = zeros(FloatType, Ni,Nj,Nk)
      z = zeros(FloatType, Ni,Nj,Nk)
  
      for i = 1:Ni, j = 1:Nj, k = 1:Nk
        x[i,j,k] = i*i/Ni/Ni
        y[i,j,k] = sqrt(j/Nj)
        z[i,j,k] = k/Nk
      end
    end

  
    # Create some scalar and vector data.
    p = zeros(FloatType, Ni, Nj, Nk)
    q = zeros(FloatType, Ni, Nj, Nk)
    vec = zeros(FloatType, 3, Ni, Nj, Nk)
    vs = zeros(SVector{3, FloatType}, Ni, Nj, Nk)  # this is an alternative way of specifying a vector dataset

    # 3×3 tensors
    tensor = zeros(FloatType, 3, 3, Ni, Nj, Nk)
    ts = zeros(SMatrix{3, 3, FloatType, 9}, Ni, Nj, Nk)

    for k = 1:Nk, j = 1:Nj, i = 1:Ni
      p[i, j, k] = i*i + k
      q[i, j, k] = k*sqrt(j)  
      vec[1, i, j, k] = i
      vec[2, i, j, k] = j
      vec[3, i, j, k] = k
      vs[i, j, k] = (i, j, k)  

      A = similar(eltype(ts))
      c = i - 2j + 3k
      for I ∈ CartesianIndices(A)
        v = (10 + I[1] - I[2]) * c
        A[I] = v
        tensor[I, i, j, k] = v
      end
      ts[i, j, k] = A
    end

    # Create some scalar data at grid cells.
    # Note that in structured grids, the cells are the hexahedra (3D) or quads (2D)
    # formed between grid points.
    local cdata
    if dim == 2
      cdata = zeros(FloatType, Ni-1, Nj-1)
      for j = 1:Nj-1, i = 1:Ni-1
        cdata[i, j] = 2i + 20 * sin(3*pi * (j-1) / (Nj-2))
      end
    elseif dim == 3
      cdata = zeros(FloatType, Ni-1, Nj-1, Nk-1)
      for k = 1:Nk-1, j = 1:Nj-1, i = 1:Ni-1
        cdata[i, j, k] = 2i + 3k * sin(3*pi * (j-1) / (Nj-2))
      end
    end

    # Test extents (this is optional!!)
    ext = map(N -> (1:N) .+ 42, (Ni, Nj, Nk))

    # Initialise new vtr file (rectilinear grid).
    local vtk
    if dim == 2
      vtk = vtk_grid(vtk_filename_struct_noext*"_$(dim)D", x, y; extent=ext, compress = compress)
    elseif dim == 3
      vtk = vtk_grid(vtk_filename_struct_noext*"_$(dim)D", x, y, z; extent=ext, compress = compress)
    end

    # Add data.
    vtk["p_values"] = p
    vtk["q_values"] = q
    
    # Test passing the second optional argument.
    @test_throws DimensionMismatch WriteVTK.num_components(
        vec, vtk, VTKCellData())
    vtk["myVector", VTKPointData()] = vec
    vtk["mySVector", VTKPointData()] = vs
    vtk["tensor"] = tensor
    vtk["tensor.SMatrix"] = ts
    vtk["myCellData"] = cdata

    # Save and close vtk file.
    append!(outfiles, vtk_save(vtk))

    name = vtk_save(vtk)[1]

    # read the file back     
    @testset "$name compress=$compress" begin
      vtk_read = VTKFile(name)
      @testset "coordinates" begin 
        # read coordinates
        x_read, y_read, z_read = get_coordinates(vtk_read)
        
        @test sum(abs.(x-x_read)) == 0.0
        @test sum(abs.(y-y_read)) == 0.0
        if dim==3
          @test  sum(abs.(z-z_read)) == 0.0
        end
      end

      # point data 
      @testset "point data" begin 
        point_data = get_point_data(vtk_read);
        p_read     = get_data_reshaped(point_data["p_values"])
        @test p == p_read    
        q_read     = get_data_reshaped(point_data["q_values"])
        @test q == q_read  
        
        myVector  = get_data_reshaped(point_data["myVector"])
        @test vec == myVector
      end

      @testset "cell data" begin 
        cell_data = get_cell_data(vtk_read);
        
        myCellData = get_data_reshaped(cell_data["myCellData"], cell_data=true)  
        @test sum(abs.(cdata-myCellData)) == 0  
      end
    end
  end
end
  