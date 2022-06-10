"""
    get_origin(vtk_file)

Retrieve the vector of coordinates of the origin of a uniform grid from the given [`VTKFile`](@ref) file.

See also: [`VTKFile`](@ref)
"""
function get_origin(vtk_file)
  # obtain dataset
  dataset_element = get_imagedata_dataset(vtk_file)

  whole_extent = get_whole_extent(vtk_file)
  
  # obtain the origin
  origin = parse.(Float64, split(attribute(dataset_element, "Origin", required=true), ' '))

  if (whole_extent[5:6] == [0,0]) && (origin[3] == 0.0)
    deleteat!(origin, 3)
  end

  return origin
end


"""
    get_spacing(vtk_file)

Retrieve a vector with the regular increments in each coordinate direction of the uniform grid from the given [`VTKFile`](@ref) file.

See also: [`VTKFile`](@ref)
"""
function get_spacing(vtk_file)
  # obtain dataset
  dataset_element = get_imagedata_dataset(vtk_file)

  # obtain the spacing
  spacing = parse.(Float64, split(attribute(dataset_element, "Spacing", required=true), ' '))

  if length(get_origin(vtk_file)) == 2
    deleteat!(spacing, 3)
  end

  return spacing
end


"""
    get_whole_extent(vtk_file)

Retrieve a vector with the `WholeExtent` 6-entry vector from the uniform grid [`VTKFile`](@ref) file.

See also: [`VTKFile`](@ref)
"""
function get_whole_extent(vtk_file)
  # obtain dataset
  dataset_element = get_imagedata_dataset(vtk_file)

  # obtain extent
  whole_extent = parse.(Int, split(attribute(dataset_element, "WholeExtent", required=true), ' '))

  return whole_extent
end


"""
    get_imagedata_dataset(vtk_file)

Retrieve ImageData dataset from the given [`VTKFile`](@ref) file.

See also: [`VTKFile`](@ref)
"""
function get_imagedata_dataset(vtk_file)
    # check imagedata
    if vtk_file.file_type != "ImageData"
      error("the file_type must be ImageData.")
    end
    
    # open the file and locate the ImageData section
    root = LightXML.root(vtk_file.xml_file)
    dataset_element = root["ImageData"][1]

    return dataset_element
end