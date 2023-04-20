module ReadVTK


using Base64: base64decode
using Downloads: download

using CodecZlib: ZlibDecompressor
using LightXML: LightXML, XMLDocument, XMLElement, parse_string, attribute, has_attribute,
                child_elements, free, content, find_element

using Reexport: @reexport
@reexport using VTKBase: VTKBase

export VTKFile, VTKData, VTKDataArray, VTKCells, VTKPrimitives,           # structs
       PVTKFile, PVTKData, PVTKDataArray, PVDFile,
       get_point_data, get_cell_data, get_data, get_data_reshaped,        # get data functions
       get_points, get_cells, get_origin, get_spacing, get_primitives,    # get geometry functions
       get_coordinates, get_coordinate_data,                              # get geometry functions
       get_example_file                                                   # other functions

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
  function VTKFile(filename, xml_file, file_type, version, byte_order, compressor,
                   appended_data, n_points, n_cells)
    # Create new object
    vtk_file = new(filename, xml_file, file_type, version, byte_order, compressor,
                   appended_data, n_points, n_cells)

    # Create finalizer that releases memory for VTK file
    f(vtk_file) = free(vtk_file.xml_file)

    # Register finalizer
    return finalizer(f, vtk_file)
  end
end

# Header type is hardcoded and corresponds to VTK XML version 1.0
header_type(::VTKFile) = UInt64

# Return true if data is compressed (= XML attribute `compressor` is non-empty in VTK file)
is_compressed(vtk_file::VTKFile) = !isempty(vtk_file.compressor)


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
    xml_file_contents = (raw_file_contents[1:(offset_begin - 1)] *
                         "\n  </AppendedData>\n</VTKFile>")
  end

  # Open file and ensure that it is a valid VTK file
  xml_file = parse_string(xml_file_contents)
  root = LightXML.root(xml_file)
  @assert LightXML.name(root) == "VTKFile"

  # Extract attributes (use `required=true` to fail fast & hard in case of unexpected content)
  file_type = attribute(root, "type", required = true)
  version = VersionNumber(attribute(root, "version", required = true))
  byte_order = attribute(root, "byte_order", required = true)
  header_type = attribute(root, "header_type", required = true)
  if has_attribute(root, "compressor")
    compressor = attribute(root, "compressor", required = true)
  else
    compressor = ""
  end

  # Ensure matching file types
  if !(file_type in ("UnstructuredGrid", "ImageData", "PolyData", "RectilinearGrid",
                     "StructuredGrid"))
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
    n_points = parse(Int, attribute(piece, "NumberOfPoints", required = true))
    n_cells = parse(Int, attribute(piece, "NumberOfCells", required = true))
  elseif (file_type == "ImageData" || file_type == "RectilinearGrid" ||
          file_type == "StructuredGrid")
    dataset_element = root[file_type][1]
    whole_extent = parse.(Int,
                          split(attribute(dataset_element, "WholeExtent", required = true),
                                ' '))
    n_points_per_grid_dir = [whole_extent[2 * i] + 1 for i in (1:3)]
    n_points = prod(n_points_per_grid_dir)
    n_cells = prod(n_points_per_grid_dir .- 1)
  elseif file_type == "PolyData"
    piece = root[file_type][1]["Piece"][1]
    n_points = parse(Int, attribute(piece, "NumberOfPoints", required = true))
    # TODO: decide how to handle the number of cells correctly, see
    #       https://github.com/JuliaVTK/ReadVTK.jl/pull/11
    n_cells = typemin(Int)
  end

  # Create and return VTKFile
  return VTKFile(filename, xml_file, file_type, version, byte_order, compressor,
                 appended_data, n_points, n_cells)
end

# Show basic information on REPL
function Base.show(io::IO, vtk_file::VTKFile)
  return print(io, "VTKFile(",
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
function piece(vtk_file::VTKFile)
  return LightXML.root(vtk_file.xml_file)[vtk_file.file_type][1]["Piece"][1]
end

"""
    isstructured(xml_file)

Returns `true` if it is a structured grid.
"""
function isstructured(xml_file)
  root = LightXML.root(xml_file)
  type = attribute(root, "type", required = true)
  if (type == "RectilinearGrid" || type == "ImageData" ||
      type == "PRectilinearGrid" || type == "PImageData" ||
      type == "StructuredGrid" || type == "PStructuredGrid")
    structured = true
  else
    structured = false
  end

  return structured
end


"""
    PVTKFile

Hold all relevant information about a Parallel VTK XML file that has been read in.

# Fields
- `filename`: original path to the PVTK file that has been read in
- `xml_file`: xml info
- `file_type`: currently only `"PRectilinearGrid"` or `"PImageData"` are supported
- `version`: currently only v1.0 is supported
- `vtk_filenames`: vector with strings that contain the filenames of each of the parallel files
- `vtk_files`: vector with `VTKFile` data that contains the info about each of the files
"""
struct PVTKFile
  filename::String
  xml_file::XMLDocument
  file_type::String
  version::VersionNumber
  vtk_filenames::Vector{String}
  vtk_files::Vector{VTKFile}
end

"""
    PVTKFile(filename; dir="")

Read in and parse the PVTK XML file specified by its `filename`.
Optionally, an additional directory name `dir` can be specified for the
location of the underlying (serial) VTK files.
"""
function PVTKFile(filename; dir = "")
  # Read in file into memory as a string
  xml_file_contents = read(filename, String)

  # Check if file begins with string that indicates this is a VTK file but *not* in XML format
  if startswith(xml_file_contents, "# pvtk DataFile Version")
    error("bad PVTK file format (found legacy format; only PVTK XML files are currently supported)")
  end

  # Open file and ensure that it is a valid VTK file
  xml_file = LightXML.parse_string(xml_file_contents)
  root = LightXML.root(xml_file)
  @assert LightXML.name(root) == "VTKFile"

  # Extract attributes (use `required=true` to fail fast & hard in case of unexpected content)
  file_type = attribute(root, "type", required = true)
  version = VersionNumber(attribute(root, "version", required = true))

  # Ensure matching file types
  if !(file_type in ("PImageData", "PRectilinearGrid", "PUnstructuredGrid",
                     "PStructuredGrid"))
    error("Unsupported file type: ", file_type)
  end

  # Ensure correct version
  @assert version == v"1.0"

  # Extract names of files & load the data
  pieces = root[file_type][1]["Piece"]
  n_pieces = length(pieces)
  vtk_filenames = Vector{String}(undef, n_pieces)
  vtk_files = Vector{VTKFile}(undef, n_pieces)

  for i in 1:n_pieces
    vtk_filenames[i] = attribute(pieces[i], "Source", required = true)
    vtk_files[i] = VTKFile(joinpath(dir, vtk_filenames[i]))
  end

  return PVTKFile(filename, xml_file, file_type, version, vtk_filenames, vtk_files)
end

# Reduce noise:
function Base.show(io::IO, vtk_file::PVTKFile)
  return print(io, "PVTKFile()")
end

"""
    PVDFile

Hold all relevant information about a PVD file that has been read in.

# Fields
- `filename`: original path to the PVTK file that has been read in
- `file_type`: currently only `"PRectilinearGrid"` or `"PImageData"` are supported
- `vtk_filenames`: vector with strings that contain the filenames of each of the files
- `directories`: vector with strings that contain the directories where each of the files are
- `timestep`: vector with `Float64` that contains the time of each of the files
"""
struct PVDFile
  filename::String
  file_type::String
  vtk_filenames::Vector{String}
  directories::Vector{String}
  timesteps::Vector{Float64}
end

"""
    PVDFile(filename)

Read in and parse the PVD XML file specified by its `filename`.
"""
function PVDFile(filename)
  # Read in file into memory as a string
  xml_file_contents = read(filename, String)

  # Open file and ensure that it is a valid VTK file
  xml_file = LightXML.parse_string(xml_file_contents)
  root = LightXML.root(xml_file)
  @assert LightXML.name(root) == "VTKFile"

  # Extract attributes (use `required=true` to fail fast & hard in case of unexpected content)
  file_type = attribute(root, "type", required = true)
  version = VersionNumber(attribute(root, "version", required = true))

  # Ensure matching file types
  if file_type != "Collection"
    error("Unsupported PVD file type: ", file_type)
  end

  # Ensure correct version
  @assert version == v"1.0"

  # Extract names of files & load the data
  pieces = root[file_type][1]["DataSet"]
  n_pieces = length(pieces)
  vtk_filenames = Vector{String}(undef, n_pieces)
  directories = Vector{String}(undef, n_pieces)
  timesteps = Vector{Float64}(undef, n_pieces)

  for i in 1:n_pieces
    file_dir = attribute(pieces[i], "file", required = true)
    vtk_filenames[i] = file_dir
    directories[i] = dirname(file_dir)
    timesteps[i] = parse(Float64, attribute(pieces[i], "timestep", required = true))
  end

  return PVDFile(filename, file_type, vtk_filenames, directories, timesteps)
end

# Reduce noise:
function Base.show(io::IO, d::PVDFile)
  return print(io, "PVDFile()")
end


include("get_functions.jl")

# Auxiliary methods
Base.keys(data::PVTKFile) = tuple(data.vtk_filenames...)
Base.length(pvtk_file::PVTKFile) = length(pvtk_file.vtk_filenames)

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

"""
    PVTKData

Convenience type to hold information about available data arrays for cells or points.

Supports a collectible/dictionary-like syntax, e.g., `keys(pvtk_data)` to show available data arrays
or `pvtk_data["varname"]` to retrieve the `VTKDataArray` for variable `varname`.
"""
struct PVTKData
  parent_xml::XMLDocument
  data::Vector{VTKData}
end

# Reduce REPL noise by defining `show`
Base.show(io::IO, data::PVTKData) = print(io, "PVTKData()")


# Retrieve a data section (should be `CellData` or `PointData`) from the VTK file
function get_data_section(vtk_file::VTKFile, section)
  names = String[]
  data_arrays = XMLElement[]

  # Iterate over XML elemens in the section
  for xml_element in child_elements(piece(vtk_file)[section][1])
    # We do not know how to handle anything other than `DataArray`s
    @assert LightXML.name(xml_element) == "DataArray"

    # Store the name and the XML element for each found data array
    push!(names, attribute(xml_element, "Name", required = true))
    push!(data_arrays, xml_element)
  end

  return VTKData(names, data_arrays, vtk_file)
end


"""
    get_cell_data(vtk_file::VTKFile)

Retrieve a lightweight `VTKData` object with the cell data of the given VTK file.

See also: [`VTKData`](@ref), [`get_point_data`](@ref)
"""
get_cell_data(vtk_file::VTKFile) = get_data_section(vtk_file, "CellData")

"""
    get_cell_data(pvtk_file::PVTKFile)

Retrieve a lightweight vector with `PVTKData` objects with the cell data of the given PVTK files.

See also: [`PVTKData`](@ref), [`get_cell_data`](@ref)
"""
function get_cell_data(pvtk_file::PVTKFile)
  n_files = length(pvtk_file)
  cdata_v = Vector{VTKData}(undef, n_files)
  for i in 1:n_files
    cdata_v[i] = get_cell_data(pvtk_file.vtk_files[i])
  end

  return PVTKData(pvtk_file.xml_file, cdata_v)
end


"""
    get_point_data(vtk_file::VTKFile)

Retrieve a lightweight `VTKData` object with the point data of the given VTK file.

See also: [`VTKData`](@ref), [`get_cell_data`](@ref)
"""
get_point_data(vtk_file::VTKFile) = get_data_section(vtk_file, "PointData")

"""
    get_point_data(pvtk_file::PVTKFile)

Retrieve a lightweight vector with `PVTKData` objects with the point data of the given PVTK files.

See also: [`PVTKData`](@ref), [`get_cell_data`](@ref)
"""
function get_point_data(pvtk_file::PVTKFile)
  n_files = length(pvtk_file)
  pdata_v = Vector{VTKData}(undef, n_files)
  for i in 1:n_files
    pdata_v[i] = get_point_data(pvtk_file.vtk_files[i])
  end

  return PVTKData(pvtk_file.xml_file, pdata_v)
end

"""
    get_coordinate_data(vtk_file::VTKFile)

Retrieve a lightweight `VTKData` object with the coordinate data of the given VTK file.

See also: [`VTKData`](@ref), [`get_point_data`](@ref),  [`get_cell_data`](@ref)
"""
get_coordinate_data(vtk_file::VTKFile) = get_data_section(vtk_file, "Coordinates")

"""
    get_coordinate_data(pvtk_file::PVTKFile)

Retrieve a lightweight `{VTKData` object with the coordinate data of the given VTK file.

See also: [`PVTKData`](@ref), [`get_point_data`](@ref),  [`get_cell_data`](@ref)
"""
function get_coordinate_data(pvtk_file::PVTKFile)
  return get_data_section.(pvtk_file.vtk_files, "Coordinates")
end

# Auxiliary methods for conveniently using `VTKData` objects like a dictionary/collectible
Base.firstindex(data::VTKData) = first(data.names)
Base.lastindex(data::VTKData) = last(data.names)
Base.length(data::VTKData) = length(data.names)
Base.size(data::VTKData) = (length(data),)
Base.keys(data::VTKData) = tuple(data.names...)

Base.firstindex(data::PVTKData) = first(data.data[1].names)
Base.lastindex(data::PVTKData) = last(data.data[1].names)
Base.length(data::PVTKData) = length(data.data[1].names)
Base.size(data::PVTKData) = (length(data.data[1]),)
Base.keys(data::PVTKData) = tuple(data.data[1].names...)

function Base.iterate(data::VTKData, state = 1)
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


function Base.getindex(data::PVTKData, name)
  index = findfirst(isequal(name), data.data[1].names)

  if isnothing(index)
    throw(KeyError(name))
  end

  n_datasets = length(data.data)
  data_array = Vector{VTKDataArray}(undef, n_datasets)
  for i in 1:n_datasets
    data_array[i] = VTKDataArray(data.data[i].data_arrays[index], data.data[i].vtk_file)
  end

  return PVTKDataArray(data.parent_xml, data_array)
end

# Reduce noise:
Base.eltype(::VTKData) = Pair{String, VTKDataArray}


"""
    VTKDataArray{T, N, Format}

Hold information about a single VTK data array (cell data or point data).

The data type is encoded as `T`, `N` represents the size of the second dimension in case of
multi-dimensonal arrays (or `1` for a vector), and `Format` encodes in which format the data is
stored in the XML file.

The actual data can be retrieved by calling `get_data` on the `PVTKDataArray` object.

See also: [`get_data`](@ref)
"""
struct VTKDataArray{T, N, Format}
  name::String
  offset::Int
  data_array::XMLElement
  vtk_file::VTKFile
end

# convencience functions
number_of_components(d::VTKDataArray{T, N, Format}) where {T, N, Format} = N
datatype(d::VTKDataArray{T, N, Format}) where {T, N, Format} = T


"""
    PVTKDataArray

Hold information about a parallel PVTK data array (cell data or point data).
The actual data can be retrieved by calling `get_data` on the `VTKDataArray` object.

See also: [`get_data`](@ref)
"""
struct PVTKDataArray
  parent_xml::XMLDocument
  data::Vector{VTKDataArray}
end

Base.show(io::IO, vtk_file::PVTKDataArray) = print(io, "PVTKDataArray()")

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
function VTKDataArray(xml_element, vtk_file::VTKFile)
  # Ensure the correct type of the XML element
  @assert LightXML.name(xml_element) == "DataArray"

  # Extract information about the underlying data
  data_type = string_to_data_type(attribute(xml_element, "type", required = true))
  name = attribute(xml_element, "Name", required = true)
  n_components = parse(Int, attribute(xml_element, "NumberOfComponents", required = true))
  format_string = attribute(xml_element, "format", required = true)

  # An offset is only used when the format is `appended`
  if has_attribute(xml_element, "offset")
    @assert format_string == "appended"
    offset = parse(Int, attribute(xml_element, "offset", required = true))
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

  return VTKDataArray{data_type, n_components, format}(name, offset, xml_element, vtk_file)
end

# Reduce REPL noise by defining `show`
function Base.show(io::IO, data_array::VTKDataArray)
  return print(io, "VTKDataArray(\"", data_array.name, "\")")
end

# Return true if data is compressed (= XML attribute `compressor` is non-empty in VTK file)
is_compressed(data_array::VTKDataArray) = is_compressed(data_array.vtk_file)


"""
    get_data(data_array::VTKDataArray)

Retrieve actual data from a `VTKDataArray` as a one- or two-dimensional array-like container.

Note: This function is not type stable but could be - help wanted!
"""
function get_data end

# Retrieve actual data for XML data array (version for storage format "binary")
function get_data(data_array::VTKDataArray{T, N, <:FormatBinary}) where {T, N}
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
function get_data(data_array::VTKDataArray{T, N, <:FormatAppended}) where {T, N}
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

    # Pass data through ZLib decompressor
    if last > length(raw)
      @show data_array data_array.vtk_file.xml_file
      @show first, last, size(raw)
      error("mistake in get_data")
    end

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
    get_data(data_array::PVTKDataArray)

Retrieve actual data from a `PVTKDataArray` as a one- or two-dimensional array-like container.
"""
function get_data(data_array::PVTKDataArray)
  n_datasets = length(data_array.data)
  dat = Vector{Array}(undef, n_datasets)
  for i in 1:n_datasets
    dat[i] = get_data(data_array.data[i])
  end

  return dat
end

"""
    get_data_reshaped(data_array::VTKDataArray; cell_data=false)

Retrieve actual data from a `VTKDataArray` and reshape it as 1D, 2D, or 3D arrays, in case we deal with structured grids.
Note that vectors or tensors will have their components stored in the first dimension of the array.
As there is no way to automatically tell from the VTK file format whether it is a tensor, the user has to reshape this manually.
"""
function get_data_reshaped(data_array::VTKDataArray{T, N}; cell_data = false) where {T, N}
  data = get_data(data_array)

  # Retrieve size of grid
  local_size, _ = get_wholeextent(data_array.vtk_file.xml_file, cell_data)

  # reshape
  if N == 1
    data_reshaped = reshape(data, local_size...)
  else
    data_reshaped = reshape(data, N, local_size...)
  end

  return data_reshaped
end


"""

Retrieve actual data from a `PVTKDataArray` and reshapes it as 1D, 2D, or 3D arrays, in case we deal with parallel structured grids.
It also puts it in the correct location the the full grid
"""
function get_data_reshaped(data_array::PVTKDataArray; cell_data = false)
  wholeextent, min_extent = get_wholeextent(data_array.parent_xml) # global grid size
  extents = get_extents(data_array.parent_xml, min_extent)         # local extents

  if cell_data
    wholeextent = wholeextent .- 1
  end
  N = number_of_components(data_array.data[1])
  type = datatype(data_array.data[1])

  # initialize full grid
  if N == 1
    data = zeros(type, wholeextent...)
  else
    data = zeros(type, N, wholeextent...)
  end

  # collect parts
  for i in 1:length(data_array.data)
    ex = extents[i]
    if cell_data
      ex = (ex[1][1:(end - 1)], ex[2][1:(end - 1)], ex[3][1:(end - 1)])
    end
    if N == 1
      data[ex...] .= get_data_reshaped(data_array.data[i], cell_data = cell_data)
    else
      data[1:N, ex...] .= get_data_reshaped(data_array.data[i], cell_data = cell_data)
    end
  end

  return data
end


"""
  get_wholeextent(xml_file, cell_data=false)

Retrieve the size of a structured grid (ImageData, RectilinearGrid). Note that this always returns three dimensions, even if the data is 1D or 2D.
"""
function get_wholeextent(xml_file, cell_data = false)

  if !isstructured(xml_file)
    error("Only works for structured grids ")
  end

  root = LightXML.root(xml_file)
  file_type = attribute(root, "type", required = true)
  dataset_element = root[file_type][1]
  whole_extent = parse.(Int,
                        split(attribute(dataset_element, "WholeExtent", required = true),
                              ' '))

  min_ex = whole_extent[1:2:(end - 1)]
  max_ex = whole_extent[2:2:end]

  local_size = [length(min_ex[i]:max_ex[i]) for i in 1:3]

  if cell_data
    local_size = local_size .- 1
    local_size[local_size .== 0] .= 1
  end

  return local_size, min_ex
end

"""
    get_extents(xml_file, min_extent=[0;0;0])

Retrieve the local size of pieces of a structured grid (ImageData, RectilinearGrid). Note that this always returns three dimensions, even if the data is 1D or 2D.
"""
function get_extents(xml_file, min_extent = [0; 0; 0])
  if !isstructured(xml_file)
    error("Only works for structured grids ")
  end

  # Extract names of files & load the data
  root = LightXML.root(xml_file)
  file_type = attribute(root, "type", required = true)
  pieces = root[file_type][1]["Piece"]
  n_pieces = length(pieces)

  # Retrieve number of points
  extents = Vector{NTuple{3, UnitRange{Int64}}}(undef, n_pieces)
  for i in 1:n_pieces
    ex = parse.(Int, split(attribute(pieces[i], "Extent", required = true)))

    # julia starts @ 1; sometimes the minimum extent starts @ zero and sometimes @ a custom value
    ex[1:2] = ex[1:2] .- min_extent[1] .+ 1
    ex[3:4] = ex[3:4] .- min_extent[2] .+ 1
    ex[5:6] = ex[5:6] .- min_extent[3] .+ 1

    extents[i] = (ex[1]:ex[2], ex[3]:ex[4], ex[5]:ex[6])
  end

  return extents
end



"""
    get_points(vtk_file::VTKFile)

Retrieve VTK points as a two-dimensional array-like container.

The points are stored in dimension × points format. Note that in VTK, points are always stored
three-dimensional, even for 1D or 2D files.

See also: [`get_cells`](@ref)
"""
function get_points(vtk_file::VTKFile)
  # Retrieve `Points` section
  points = find_element(piece(vtk_file), "Points")
  @assert !isnothing(points)

  # Get the first `DataArray` XML element (there should be only one)
  xml_data_array = find_element(points, "DataArray")
  @assert !isnothing(xml_data_array)

  # Create data array and return actual data
  data_array = VTKDataArray(xml_data_array, vtk_file)

  return get_data(data_array)
end

"""
  get_points(vtk_file::PVTKFile)

Retrieve VTK points as a two-dimensional array-like container for a parallel file

"""
get_points(pvtk_file::PVTKFile) = get_points.(pvtk_file.vtk_files)

"""
    get_coordinates(vtk_file::VTKFile; x_string="x", y_string="y", z_string="z")

Retrieve VTK coordinate vectors in each direction as a tuple of 1D vectors for a RectilinearGrid file.

Note that in VTK, points are always stored three-dimensional, even for 1D or 2D files, so you will always retrieve a tuple with three vectors.

See also: [`get_cells`](@ref)
"""
function get_coordinates(vtk_file::VTKFile; x_string = "x", y_string = "y", z_string = "z")
  if vtk_file.file_type == "RectilinearGrid"
    coordinates = get_coordinate_data(vtk_file)
    x = get_data(coordinates[x_string])
    y = get_data(coordinates[y_string])
    z = get_data(coordinates[z_string])
  elseif vtk_file.file_type == "StructuredGrid"
    wholeextent, _ = get_wholeextent(vtk_file.xml_file)
    points = get_points(vtk_file)
    data = reshape(points, 3, wholeextent...)
    x = data[1, :, :, :]
    y = data[2, :, :, :]
    z = data[3, :, :, :]
  else
    error("The file type of the VTK file must be 'RectilinearGrid' or 'StructuredGrid' (current: $(vtk_file.file_type)).")
  end

  return x, y, z
end


"""
    get_coordinates(pvtk_file::{VTKFile; x_string="x", y_string="y", z_string="z")

Retrieve VTK coordinate vectors in each direction as a tuple of 1D vectors for a PRectilinearGrid file.

Note that in VTK, points are always stored three-dimensional, even for 1D or 2D files, so you will always retrieve a tuple with three vectors.

See also: [`get_cells`](@ref)
"""
function get_coordinates(pvtk_file::PVTKFile; x_string = "x", y_string = "y",
                         z_string = "z")
  if pvtk_file.file_type == "PStructuredGrid"
    points = get_points(pvtk_file)

    wholeextent, min_extent = get_wholeextent(pvtk_file.xml_file)  # global grid size
    extents = get_extents(pvtk_file.xml_file, min_extent)         # local extents

    type = typeof(points[1][1])

    # initialize full grid
    data = zeros(type, 3, wholeextent...)

    # collect parts
    for i in 1:length(points)
      ex = extents[i]
      data[1:3, ex...] .= reshape(points[i], 3, length.(extents[i])...)
      x = data[1, :, :, :]
      y = data[2, :, :, :]
      z = data[3, :, :, :]
    end
  elseif pvtk_file.file_type == "PRectilinearGrid"
    coords = get_coordinates.(pvtk_file.vtk_files, x_string = x_string, y_string = y_string,
                              z_string = z_string)
    x, y, z = coords[1][1][:], coords[1][2][:], coords[1][3][:]
    for i in 2:length(pvtk_file)
      x = [x; coords[i][1][:]]
      y = [y; coords[i][2][:]]
      z = [z; coords[i][3][:]]
    end
    x, y, z = unique(x), unique(y), unique(z)

  else
    error("File should be of type PRectilinearGrid or PStructuredGrid. Current file type: $(pvtk_file.file_type)")
  end

  return x, y, z
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
function get_cells(vtk_file::VTKFile)
  # Retrieve `Cells` section
  cells = find_element(piece(vtk_file), "Cells")
  @assert !isnothing(cells)

  # Iterate over available XML data arrays and store corresponding data in `VTKDataArray`s
  connectivity = offsets = types = nothing
  for xml_element in cells["DataArray"]
    a = attribute(xml_element, "Name", required = true)
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

  # Create VTKCells container and convert VTK's zero-based indices to Julia's one-based indices
  # Note: the offsets do not need to be updated since they point *past* the last entry
  # in VTK files (C-style), while in Julia it is custom to point *at* the last entry.
  return VTKCells(get_data(connectivity) + oneunit.(get_data(connectivity)),
                  get_data(offsets),
                  get_data(types))
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
function get_primitives(vtk_file::VTKFile, primitive_type::AbstractString)
  @assert vtk_file.file_type == "PolyData"
  if !(primitive_type in ("Verts", "Lines", "Polys"))
    error("Unsupported `primitive type`: \"", primitive_type,
          "\". Supported values are: \"Verts\", \"Lines\", or \"Polys\".")
  end

  xml = find_element(piece(vtk_file), primitive_type)
  @assert !isnothing(xml) string("This PolyData file has no primitives of type \"",
                                 primitive_type, "\".")

  # Iterate over available XML data arrays and store corresponding data in `VTKDataArray`s
  connectivity = nothing
  offsets = nothing
  for xml_element in xml["DataArray"]
    a = attribute(xml_element, "Name", required = true)
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
  return VTKPrimitives(get_data(connectivity), get_data(offsets))
end

# Convenience functions for working with `VTKPrimitives` container
Base.length(primitives::VTKPrimitives) = length(primitives.offsets)
Base.size(primitives::VTKPrimitives) = (length(primitives),)


"""
    get_example_file(filename; head="main", output_directory=".", force=false)

Retrieve an example file from the
[`ReadVTK_examples` repository](https://github.com/JuliaVTK/ReadVTK_examples)
at commit/branch `head` and store it in the `output_directory`. If the file already
exists locally, do not download the file again unless `force` is true. Return the local path to the
downloaded file.
"""
function get_example_file(filename; head = "main", output_directory = ".", force = false)
  filepath = joinpath(output_directory, filename)
  if !isfile(filepath) || force
    url = ("https://github.com/JuliaVTK/ReadVTK_examples/raw/"
           * head
           * "/examples/"
           * filename)
    download(url, filepath)
  end

  return filepath
end

end # module
