"""
    get_origin(vtk_file)

Retrieve the vector of coordinates of the origin of a uniform grid from the given [`VTKFile`](@ref) file.

See also: [`VTKFile`](@ref)
"""
function get_origin(vtk_file)
  # check imagedata
  if vtk_file.file_type != "ImageData"
    error("origin can be read only for ImageData file_type.")
  end

  # open the file and locate the ImageData section
  root = LightXML.root(vtk_file.xml_file)
  dataset_element = root["ImageData"][1]

  # obtain the origin
  origin = parse.(Float64, split(attribute(dataset_element, "Origin", required=true), ' '))

  return origin
end



"""
    get_spacing(vtk_file)

Retrieve a vector with the regular increments in each coordinate direction of the uniform grid from the given [`VTKFile`](@ref) file.

See also: [`VTKFile`](@ref)
"""
function get_spacing(vtk_file)
  # check imagedata
  if vtk_file.file_type != "ImageData"
    error("spacing can be read only for ImageData file_type.")
  end
  
  # open the file and locate the ImageData section
  root = LightXML.root(vtk_file.xml_file)
  dataset_element = root["ImageData"][1]

  # obtain the origin
  spacing = parse.(Float64, split(attribute(dataset_element, "Spacing", required=true), ' '))

  return spacing
end

