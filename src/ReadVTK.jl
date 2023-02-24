module ReadVTK


using Base64: base64decode
using Downloads: download

using CodecZlib: ZlibDecompressor
using LightXML: LightXML, XMLDocument, XMLElement, parse_string, attribute, has_attribute,
                child_elements, free, content, find_element

export VTKFile, VTKData, VTKDataArray, VTKCells, VTKPrimitives,        # structs
       get_point_data, get_cell_data, get_data, get_data_reshaped,     # get data functions
       get_points, get_cells, get_origin, get_spacing, get_primitives, # get geometry functions
       get_coordinates, get_coordinate_data,                           # get geometry functions
       get_example_file                                                # other functions

"""
    VTKFile

Hold all relevant information about a VTK XML file that has been read in.

# Fields
- `filename`: original path to the VTK file that has been read in
- `xml_file`: object that represents the XML file
- `file_type`: currently only `"UnstructuredGrid"` or `"ImageData"` are supported
- `version`: currently only v1.0 is supported
- `byte_order`: can be `LittleEndian` or `BigEndian` and must currently be the same as the system's
- `compressor`: can be empty (no compression) or `vtkZLibDataCompressor`
- `appended_data`: in case of appended data (see XML documentation), the data is stored here for
                   convenient retrieval (otherwise it is empty)
- `n_points`: number of points in the VTK file
- `n_cells`: number of cells in the VTK file`
"""
mutable struct VTKFile
  filename::String
  xml_file::XMLDocument
  file_type::String
  version::VersionNumber
  byte_order::String
  compressor::String
  appended_data::Vector{UInt8}
  n_points::Int
  n_cells::Int

  # Create inner constructor to add finalizer (this requires VTKFile to be a *mutable* struct)
  function VTKFile(filename, xml_file, file_type, version, byte_order, compressor, appended_data,
                   n_points, n_cells)
    # Create new object
    vtk_file = new(filename, xml_file, file_type, version, byte_order, compressor, appended_data,
                   n_points, n_cells)

    # Create finalizer that releases memory for VTK file
    f(vtk_file) = free(vtk_file.xml_file)

    # Register finalizer
    finalizer(f, vtk_file)
  end
end

# Header type is hardcoded and corresponds to VTK XML version 1.0
header_type(::VTKFile) = UInt64


# Return true if data is compressed (= XML attribute `compressor` is non-empty in VTK file)
is_compressed(vtk_file::VTKFile) = !isempty(vtk_file.compressor)

include("get_functions.jl")

"""
    VTKFile(filename)

Read in and parse the VTK XML file specified by its `filename`.
"""
function VTKFile(filename)
  # Read in file into memory as a string
  raw_file_contents = read(filename, String)

  # Check if file begins with string that indicates this is a VTK file but *not* in XML format
  if startswith(raw_file_contents, "# vtk DataFile Version")
    error("bad VTK file format (found legacy format; only VTK XML files are currently supported)")
  end

  # Check if data is stored in "appended" mode
  marker = findfirst("<AppendedData encoding=\"raw\">", raw_file_contents)
  if isnothing(marker)
    # `appended_data` remains empty and the LightXML library may parse the entire string
    appended_data = UInt8[]
    xml_file_contents = raw_file_contents
  else
    # Since data is stored in appended mode, the VTK file is technically not valid XML. Therefore we
    # extract the appended data and store it as a byte stream in `appended_data`. We then clean up
    # the string (i.e., remove appended data) and pass a valid XML file string to the LightXML
    # library
    offset_begin = first(findnext("_", raw_file_contents, last(marker))) + 1
    offset_end = first(findnext("</AppendedData>", raw_file_contents, offset_begin)) - 1
    appended_data = Vector{UInt8}(rstrip(raw_file_contents[offset_begin:offset_end]))
    xml_file_contents = raw_file_contents[1:offset_begin-1] * "\n  </AppendedData>\n</VTKFile>"
  end

  # Open file and ensure that it is a valid VTK file
  xml_file = parse_string(xml_file_contents)
  root = LightXML.root(xml_file)
  @assert LightXML.name(root) == "VTKFile"

  # Extract attributes (use `required=true` to fail fast & hard in case of unexpected content)
  file_type = attribute(root, "type", required=true)
  version = VersionNumber(attribute(root, "version", required=true))
  byte_order = attribute(root, "byte_order", required=true)
  header_type = attribute(root, "header_type", required=true)
  if has_attribute(root, "compressor")
    compressor = attribute(root, "compressor", required=true)
  else
    compressor = ""
  end

  # Ensure matching file types
  if !(file_type in ("UnstructuredGrid", "ImageData", "PolyData", "RectilinearGrid"))
    error("Unsupported file type: ", file_type)
  end

  # Ensure correct version
  @assert version == v"1.0"

  # Ensure matching byte order
  is_little_endian = ENDIAN_BOM == 0x04030201
  @assert ((byte_order == "LittleEndian" && is_little_endian) ||
           (byte_order == "BigEndian" && !is_little_endian))

  # Ensure matching header type
  @assert header_type == "UInt64"

  # Ensure supported compression type
  @assert in(compressor, ("", "vtkZLibDataCompressor"))

  if file_type == "UnstructuredGrid"
    # Extract piece
    piece = root[file_type][1]["Piece"][1]
    n_points = parse(Int, attribute(piece, "NumberOfPoints", required=true))
    n_cells = parse(Int, attribute(piece, "NumberOfCells", required=true))
  elseif file_type == "ImageData" || file_type == "RectilinearGrid"
    dataset_element = root[file_type][1]
    whole_extent = parse.(Int, split(attribute(dataset_element, "WholeExtent", required=true), ' '))
    n_points_per_grid_dir = [whole_extent[2*i]+1 for i in (1:3)]
    n_points = prod(n_points_per_grid_dir)
    n_cells = prod(n_points_per_grid_dir .- 1)
  elseif file_type == "PolyData"
    piece = root[file_type][1]["Piece"][1]
    n_points = parse(Int, attribute(piece, "NumberOfPoints", required=true))
    # TODO: decide how to handle the number of cells correctly, see  
    #       https://github.com/trixi-framework/ReadVTK.jl/pull/11
    n_cells = typemin(Int)
  end

  # Create and return VTKFile
  VTKFile(filename, xml_file, file_type, version, byte_order, compressor, appended_data, n_points,
          n_cells)
end

# Show basic information on REPL
function Base.show(io::IO, vtk_file::VTKFile)
  print(io, "VTKFile(",
            "\"", vtk_file.filename, "\", ",
            "<XMLDocument>", ", ",
            "\"", vtk_file.file_type, "\", ",
            "\"", vtk_file.version, "\", ",
            "\"", vtk_file.byte_order, "\", ",
            "\"", vtk_file.compressor, "\", ",
            "<appended_data>", ", ",
            vtk_file.n_points, ", ",
            vtk_file.n_cells, ")")
end

# Return `Piece` XML element that contains all VTK data
piece(vtk_file::VTKFile) = LightXML.root(vtk_file.xml_file)[vtk_file.file_type][1]["Piece"][1]


"""
    VTKData

Convenience type to hold information about available data arrays for cells or points.

Supports a collectible/dictionary-like syntax, e.g., `keys(vtk_data)` to show available data arrays
or `vtk_data["varname"]` to retrieve the `VTKDataArray` for variable `varname`.
"""
struct VTKData
  names::Vector{String}
  data_arrays::Vector{XMLElement}
  vtk_file::VTKFile
end

# Reduce REPL noise by defining `show`
Base.show(io::IO, data::VTKData) = print(io, "VTKData()")

# Retrieve a data section (should be `CellData` or `PointData`) from the VTK file
function get_data_section(vtk_file, section)
  names = String[]
  data_arrays = XMLElement[]

  # Iterate over XML elemens in the section
  for xml_element in child_elements(piece(vtk_file)[section][1])
    # We do not know how to handle anything other than `DataArray`s
    @assert LightXML.name(xml_element) == "DataArray"

    # Store the name and the XML element for each found data array
    push!(names, attribute(xml_element, "Name", required=true))
    push!(data_arrays, xml_element)
  end

  VTKData(names, data_arrays, vtk_file)
end


"""
    get_cell_data(vtk_file)

Retrieve a lightweight `VTKData` object with the cell data of the given VTK file.

See also: [`VTKData`](@ref), [`get_point_data`](@ref)
"""
get_cell_data(vtk_file) = get_data_section(vtk_file, "CellData")


"""
    get_point_data(vtk_file)

Retrieve a lightweight `VTKData` object with the point data of the given VTK file.

See also: [`VTKData`](@ref), [`get_cell_data`](@ref)
"""
get_point_data(vtk_file) = get_data_section(vtk_file, "PointData")

"""
    get_coordinate_data(vtk_file)

Retrieve a lightweight `VTKData` object with the coordinate data of the given VTK file.

See also: [`VTKData`](@ref), [`get_point_data`](@ref),  [`get_cell_data`](@ref)
"""
get_coordinate_data(vtk_file) = get_data_section(vtk_file, "Coordinates")


# Auxiliary methods for conveniently using `VTKData` objects like a dictionary/collectible
Base.firstindex(data::VTKData) = first(data.names)
Base.lastindex(data::VTKData) = last(data.names)
Base.length(data::VTKData) = length(data.names)
Base.size(data::VTKData) = (length(data),)
Base.keys(data::VTKData) = tuple(data.names...)

function Base.iterate(data::VTKData, state=1)
  if state > length(data)
    return nothing
  else
    return (data.names[state] => data[data.names[state]], state + 1)
  end
end

function Base.getindex(data::VTKData, name)
  index = findfirst(isequal(name), data.names)

  if isnothing(index)
    throw(KeyError(name))
  end

  return VTKDataArray(data.data_arrays[index], data.vtk_file)
end

Base.eltype(::VTKData) = Pair{String, VTKDataArray}


"""
    VTKDataArray{T, N, Format}

Hold information about a single VTK data array (cell data or point data).

The data type is encoded as `T`, `N` represents the size of the second dimension in case of
multi-dimensonal arrays (or `1` for a vector), and `Format` encodes in which format the data is
stored in the XML file.

The actual data can be retrieved by calling `get_data` on the `VTKDataArray` object.

See also: [`get_data`](@ref)
"""
struct VTKDataArray{T, N, Format}
  name::String
  offset::Int
  data_array::XMLElement
  vtk_file::VTKFile
end

# Auxiliary types for type stability
struct FormatBinary end
struct FormatAppended end
struct FormatAscii end

# Types taken from https://vtk.org/doc/release/7.0/html/vtkType_8h_source.html
function string_to_data_type(s)
  if s == "Int8"
    Int8
  elseif s == "Int16"
    Int16
  elseif s == "Int32"
    Int32
  elseif s == "Int64"
    Int64
  elseif s == "UInt8"
    UInt8
  elseif s == "UInt16"
    UInt16
  elseif s == "UInt32"
    UInt32
  elseif s == "UInt64"
    UInt64
  elseif s == "Float32"
    Float32
  elseif s == "Float64"
    Float64
  else
    error("unknown data type: ", s)
  end
end


"""
    VTKDataArray(xml_element, vtk_file)

Create a lightweight container for a given `xml_element` and the corresponding `vtk_file`.

The `xml_element` must be of type `DataArray` and the `vtk_file` parameter is required for
retrieving the actual data. You can pass a `VTKDataArray` object to `get_data` to retrieve the actual data.

See also: [`get_data`](@ref)
"""
function VTKDataArray(xml_element, vtk_file)
  # Ensure the correct type of the XML element
  @assert LightXML.name(xml_element) == "DataArray"

  # Extract information about the underlying data
  data_type = string_to_data_type(attribute(xml_element, "type", required=true))
  name = attribute(xml_element, "Name", required=true)
  n_components = parse(Int, attribute(xml_element, "NumberOfComponents", required=true))
  format_string = attribute(xml_element, "format", required=true)

  # An offset is only used when the format is `appended`
  if has_attribute(xml_element, "offset")
    @assert format_string == "appended"
    offset = parse(Int, attribute(xml_element, "offset", required=true))
  else
    offset = -1
  end

  # Convert string to type for type stability later on
  if format_string == "binary"
    format = FormatBinary
  elseif format_string == "appended"
    format = FormatAppended
  elseif format_string == "ascii"
    format = FormatAscii
  else
    error("unknown data array format: ", format_string)
  end

  VTKDataArray{data_type, n_components, format}(name, offset, xml_element, vtk_file)
end

# Reduce REPL noise by defining `show`
Base.show(io::IO, data_array::VTKDataArray) = print(io, "VTKDataArray(\"", data_array.name, "\")")

# Return true if data is compressed (= XML attribute `compressor` is non-empty in VTK file)
is_compressed(data_array::VTKDataArray) = is_compressed(data_array.vtk_file)


"""
    get_data(data_array::VTKDataArray)

Retrieve actual data from a `VTKDataArray` as a one- or two-dimensional array-like container.

Note: This function is not type stable but could be - help wanted!
"""
function get_data(data_array::VTKDataArray) end

# Retrieve actual data for XML data array (version for storage format "binary")
function get_data(data_array::VTKDataArray{T,N,<:FormatBinary}) where {T,N}
  # To retrieve the raw data,
  # * first get the content of the corresponding XML data array
  # * then remove leading/trailing whitespace
  # * finally decode from Base64 to binary representation
  raw = base64decode(strip(content(data_array.data_array)))

  if is_compressed(data_array)
    # If data is stored compressed, the first four integers of type `header_type` are the header and
    # must be discarded
    first = 1 + 4 * sizeof(header_type(data_array.vtk_file))
    last = length(raw)

    # Pass data through ZLib decompressor
    uncompressed = transcode(ZlibDecompressor, raw[first:last])
  else
    # If data is stored uncompressed, the first integer of type `header_type` is the header and must
    # be discarded
    first = 1 + sizeof(header_type(data_array.vtk_file))
    last = length(raw)
    uncompressed = view(raw, first:last)
  end

  # Convert from binary representation to actual type
  v = reinterpret(T, uncompressed)

  if N == 1
    # If it is a one dimensional array, just return it as is
    return v
  else 
    # Otherwise, reshape into two-dimensional array
    return reshape(v, N, :)
  end
end

# Retrieve actual data for XML data array (version for storage format "appended")
function get_data(data_array::VTKDataArray{T,N,<:FormatAppended}) where {T,N}
  raw = data_array.vtk_file.appended_data
  HeaderType = header_type(data_array.vtk_file)

  if is_compressed(data_array)
    # If data is stored compressed, the first four integers of type `header_type` are the header and
    # the fourth value contains the number of bytes to read
    first = data_array.offset + 1
    last = data_array.offset + 4 * sizeof(HeaderType)
    header = Int.(reinterpret(HeaderType, raw[first:last]))
    n_bytes = header[4]

    first = data_array.offset + 4 * sizeof(HeaderType) + 1
    last = first + n_bytes - 1
    uncompressed = transcode(ZlibDecompressor, raw[first:last])
  else
    # If data is stored uncompressed, the first integer of type `header_type` is the header and
    # contains the number of bytes to read
    first = data_array.offset + 1
    last = data_array.offset + sizeof(HeaderType)
    header = Int.(reinterpret(HeaderType, raw[first:last]))
    n_bytes = header[1]
    first = data_array.offset + sizeof(HeaderType) + 1
    last = first + n_bytes - 1
    uncompressed = view(raw, first:last)
  end

  v = reinterpret(T, uncompressed)
  if N == 1
    # If it is a one dimensional array, just return it as is
    return v
  else 
    # Otherwise, reshape into two-dimensional array
    return reshape(v, N, :)
  end
end


"""
    get\\_data\\_reshaped(data_array::VTKDataArray; cell_data=false)

Retrieve actual data from a `VTKDataArray` and reshapes them as 1D, 2D or 3D arrays, in case we deal with structured grids.
Note that vectors or tensors will have their components stored in the first dimension of the array. 
As there is no way to automatically tell from the VTK file format whether it is a tensor, the user has to reshape this accordingly.

"""
function get_data_reshaped(data_array::VTKDataArray{T,N}; cell_data=false) where {T,N}

  data = get_data(data_array);

  # Retrieve size of grid
  local_size = get_local_size(data_array.vtk_file.xml_file, data_array.vtk_file.file_type, cell_data)
  
  # reshape 
  if N == 1
    data_reshaped = reshape(data, local_size...)
  else
    data_reshaped = reshape(data, N, local_size...)
  end

  return data_reshaped
end

"""
    get_local_size(xml_file, file_type, cell_data=false)

Retrieve the local size of a structured grid (ImageData, RectilinearGrid). Note that this always returns three dimensions, even if the data is 1D or 2D. 
"""
function get_local_size(xml_file, file_type, cell_data=false)

  root = LightXML.root(xml_file)
  dataset_element = root[file_type][1]
  whole_extent = parse.(Int, split(attribute(dataset_element, "WholeExtent", required=true), ' '))

  N = whole_extent[2:2:end] - whole_extent[1:2:end-1] .+ 1

  if cell_data
    N = N .- 1
    N[N.==0] .= 1
  end


  return N
end

"""
    get_points(vtk_file)

Retrieve VTK points as a two-dimensional array-like container.

The points are stored in dimension × points format. Note that in VTK, points are always stored
three-dimensional, even for 1D or 2D files.

See also: [`get_cells`](@ref)
"""
function get_points(vtk_file)
  # Retrieve `Points` section
  points = find_element(piece(vtk_file), "Points")
  @assert !isnothing(points)

  # Get the first `DataArray` XML element (there should be only one)
  xml_data_array = find_element(points, "DataArray")
  @assert !isnothing(xml_data_array)

  # Create data array and return actual data
  data_array = VTKDataArray(xml_data_array, vtk_file)

  get_data(data_array)
end

"""
  get_coordinates(vtk_file; x_string="x",y_string="y",z_string="z")

Retrieve VTK coordinate vectors in each direction as 1D vectors for a RectilinearGrid file

The points are given as 1D vectors. Note that in VTK, points are always stored
three-dimensional, even for 1D or 2D files, so you will always retrieve 3 vectors.
The 

See also: [`get_cells`](@ref)
"""
function get_coordinates(vtk_file; x_string="x",y_string="y",z_string="z")

  if vtk_file.file_type != "RectilinearGrid"
      error("the file_type must be RectilinearGrid.")
  end

  coordinates = get_coordinate_data(vtk_file)
  x = get_data(coordinates[x_string])
  y = get_data(coordinates[y_string])
  z = get_data(coordinates[z_string])

  return  x,y,z
end


"""
    VTKCells{Connectivity, Offsets, Types}

Store the `connectivity`, `offsets`, and `types` information from the VTK file as one-dimensional
array-like containers. See the
[VTK file format documentation](http://vtk.org/VTK/img/file-formats.pdf) for information on how to
connect the `connectivity` and `offset` arrays to the actual geometric points.

You may use `length` to determine the number of cells from a `VTKCells` object.

See also: [`get_points`](@ref)
"""
struct VTKCells{Connectivity, Offsets, Types}
  connectivity::Connectivity
  offsets::Offsets
  types::Types
end


"""
    get_cells(vtk_file)

Retrieve VTK cell information as an object of type `VTKCells`.

See also: [`VTKCells`](@ref)
"""
function get_cells(vtk_file)
  # Retrieve `Cells` section
  cells = find_element(piece(vtk_file), "Cells")
  @assert !isnothing(cells)

  # Iterate over available XML data arrays and store corresponding data in `VTKDataArray`s
  connectivity = offsets = types = nothing
  for xml_element in cells["DataArray"]
    a = attribute(xml_element, "Name", required=true)
    if a == "connectivity"
      connectivity = VTKDataArray(xml_element, vtk_file)
    elseif a == "offsets"
      offsets = VTKDataArray(xml_element, vtk_file)
    elseif a == "types"
      types = VTKDataArray(xml_element, vtk_file)
    else
      error("unknown cell information: ", a)
    end
  end

  # Ensure that nothing is missing
  @assert !isnothing(connectivity)
  @assert !isnothing(offsets)
  @assert !isnothing(types)

  # Create VTKCells container
  VTKCells(get_data(connectivity), get_data(offsets), get_data(types))
end

# Convenience functions for working with `VTKCells` container
Base.length(vtk_cells::VTKCells) = length(vtk_cells.types)
Base.size(vtk_cells::VTKCells) = (length(vtk_cells),)


"""
    VTKPrimitives{Connectivity, Offsets}

Store the `connectivity` and `offsets` information from the VTK PolyData file as one-dimensional
array-like containers. See the
[VTK file format documentation](http://vtk.org/VTK/img/file-formats.pdf) for information on how to
connect the `connectivity` and `offset` arrays to the actual geometric points.

You may use `length` to determine the number of cells from a `VTKPrimitives` object.
"""
struct VTKPrimitives{Connectivity, Offsets}
  connectivity::Connectivity
  offsets::Offsets
end


"""
  get_primitives(vtk_file, primitive_type::AbstractString)

Retrieve VTK primitives as an object of type `VTKPrimitives`.
Supported values of `primitive type` are : \"Verts\", \"Lines\", or \"Polys\".

See also: [`VTKPrimitives`](@ref)
"""
function get_primitives(vtk_file, primitive_type::AbstractString)
  @assert vtk_file.file_type == "PolyData"
  if !(primitive_type in ("Verts", "Lines", "Polys"))
    error(
      "Unsupported `primitive type`: \"", primitive_type, "\". Supported values are: \"Verts\", \"Lines\", or \"Polys\"."
    )
  end

  xml = find_element(piece(vtk_file), primitive_type)
  @assert !isnothing(xml) string("This PolyData file has no primitives of type \"", primitive_type, "\".") 

  # Iterate over available XML data arrays and store corresponding data in `VTKDataArray`s
  connectivity = nothing
  offsets      = nothing
  for xml_element in xml["DataArray"]
    a = attribute(xml_element, "Name", required=true)
    if a == "connectivity"
      connectivity = VTKDataArray(xml_element, vtk_file)
    elseif a == "offsets"
      offsets = VTKDataArray(xml_element, vtk_file)
    end
  end

  # Ensure that nothing is missing
  @assert !isnothing(connectivity)
  @assert !isnothing(offsets)

  # Create VTKPrimitives container
  VTKPrimitives(get_data(connectivity), get_data(offsets))
end

# Convenience functions for working with `VTKPrimitives` container
Base.length(primitives::VTKPrimitives) = length(primitives.offsets)
Base.size(primitives::VTKPrimitives) = (length(primitives),)


"""
    get_example_file(filename; head="main", output_directory=".", force=false)

Retrieve an example file from the
[`ReadVTK_examples` repository](https://github.com/trixi-framework/ReadVTK_examples)
at commit/branch `head` and store it in the `output_directory`. If the file already
exists locally, do not download the file again unless `force` is true. Return the local path to the
downloaded file.
"""
function get_example_file(filename; head="main", output_directory=".", force=false)
  filepath = joinpath(output_directory, filename)
  if !isfile(filepath) || force
    url = ("https://github.com/trixi-framework/ReadVTK_examples/raw/"
           * head
           * "/examples/"
           * filename)
    download(url, filepath)
  end

  return filepath
end

end # module
