
"""
    get_origin(vtk_file)

Retrieve the origin of a VTK regular structured grid file.

"""
function get_origin(vtk_file)

  root = LightXML.root( vtk_file.xml_file )

  piece = root["ImageData"][1]

  originPoint = parse.(Float64, split( attribute(piece, "Origin",      required=true) , ' ' ) )

  return originPoint
end

