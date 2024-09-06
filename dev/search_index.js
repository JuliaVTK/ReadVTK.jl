var documenterSearchIndex = {"docs":
[{"location":"license/#License","page":"License","title":"License","text":"","category":"section"},{"location":"license/","page":"License","title":"License","text":"MIT LicenseCopyright (c) 2021-present The Trixi Authors (see https://github.com/trixi-framework/Trixi.jl)Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.","category":"page"},{"location":"contributing/#Contributing","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"ReadVTK is an open-source project and we are very happy to accept contributions from the community. Please feel free to open issues or submit patches (preferably as pull requests) any time. For planned larger contributions, it is often beneficial to get in contact with us first (e.g., by creating an issue).","category":"page"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"ReadVTK and its contributions are licensed under the MIT license (see License). As a contributor, you certify that all your contributions are in conformance with the Developer Certificate of Origin (Version 1.1), which is reproduced below.","category":"page"},{"location":"contributing/#Developer-Certificate-of-Origin-(Version-1.1)","page":"Contributing","title":"Developer Certificate of Origin (Version 1.1)","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"The following text was taken from https://developercertificate.org:","category":"page"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"Developer Certificate of Origin\nVersion 1.1\n\nCopyright (C) 2004, 2006 The Linux Foundation and its contributors.\n1 Letterman Drive\nSuite D4700\nSan Francisco, CA, 94129\n\nEveryone is permitted to copy and distribute verbatim copies of this\nlicense document, but changing it is not allowed.\n\n\nDeveloper's Certificate of Origin 1.1\n\nBy making a contribution to this project, I certify that:\n\n(a) The contribution was created in whole or in part by me and I\n    have the right to submit it under the open source license\n    indicated in the file; or\n\n(b) The contribution is based upon previous work that, to the best\n    of my knowledge, is covered under an appropriate open source\n    license and I have the right under that license to submit that\n    work with modifications, whether created in whole or in part\n    by me, under the same open source license (unless I am\n    permitted to submit under a different license), as indicated\n    in the file; or\n\n(c) The contribution was provided directly to me by some other\n    person who certified (a), (b) or (c) and I have not modified\n    it.\n\n(d) I understand and agree that this project and the contribution\n    are public and that a record of the contribution (including all\n    personal information I submit with it, including my sign-off) is\n    maintained indefinitely and may be redistributed consistent with\n    this project or the open source license(s) involved.","category":"page"},{"location":"reference/#ReadVTK.jl-API","page":"Reference","title":"ReadVTK.jl API","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"CurrentModule = ReadVTK","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [ReadVTK]","category":"page"},{"location":"reference/#ReadVTK.PVDFile","page":"Reference","title":"ReadVTK.PVDFile","text":"PVDFile\n\nHold all relevant information about a PVD file that has been read in.\n\nFields\n\nfilename: original path to the PVTK file that has been read in\nfile_type: currently only \"PRectilinearGrid\" or \"PImageData\" are supported\nvtk_filenames: vector with strings that contain the filenames of each of the files\ndirectories: vector with strings that contain the directories where each of the files are\ntimestep: vector with Float64 that contains the time of each of the files\n\n\n\n\n\n","category":"type"},{"location":"reference/#ReadVTK.PVDFile-Tuple{Any}","page":"Reference","title":"ReadVTK.PVDFile","text":"PVDFile(filename)\n\nRead in and parse the PVD XML file specified by its filename.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.PVTKData","page":"Reference","title":"ReadVTK.PVTKData","text":"PVTKData\n\nConvenience type to hold information about available data arrays for cells or points.\n\nSupports a collectible/dictionary-like syntax, e.g., keys(pvtk_data) to show available data arrays or pvtk_data[\"varname\"] to retrieve the VTKDataArray for variable varname.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ReadVTK.PVTKDataArray","page":"Reference","title":"ReadVTK.PVTKDataArray","text":"PVTKDataArray\n\nHold information about a parallel PVTK data array (cell data or point data). The actual data can be retrieved by calling get_data on the VTKDataArray object.\n\nSee also: get_data\n\n\n\n\n\n","category":"type"},{"location":"reference/#ReadVTK.PVTKFile","page":"Reference","title":"ReadVTK.PVTKFile","text":"PVTKFile\n\nHold all relevant information about a Parallel VTK XML file that has been read in.\n\nFields\n\nfilename: original path to the PVTK file that has been read in\nxml_file: xml info\nfile_type: currently only \"PRectilinearGrid\" or \"PImageData\" are supported\nversion: VTK XML file format version\nvtk_filenames: vector with strings that contain the filenames of each of the parallel files\nvtk_files: vector with VTKFile data that contains the info about each of the files\n\n\n\n\n\n","category":"type"},{"location":"reference/#ReadVTK.PVTKFile-Tuple{Any}","page":"Reference","title":"ReadVTK.PVTKFile","text":"PVTKFile(filename; dir=\"\")\n\nRead in and parse the PVTK XML file specified by its filename. Optionally, an additional directory name dir can be specified for the location of the underlying (serial) VTK files.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.VTKCells","page":"Reference","title":"ReadVTK.VTKCells","text":"VTKCells{Connectivity, Offsets, Types}\n\nStore the connectivity, offsets, and types information from the VTK file as one-dimensional array-like containers. See the VTK file format documentation for information on how to connect the connectivity and offset arrays to the actual geometric points.\n\nYou may use length to determine the number of cells from a VTKCells object.\n\nSee also: get_points\n\n\n\n\n\n","category":"type"},{"location":"reference/#ReadVTK.VTKData","page":"Reference","title":"ReadVTK.VTKData","text":"VTKData\n\nConvenience type to hold information about available data arrays for cells or points.\n\nSupports a collectible/dictionary-like syntax, e.g., keys(vtk_data) to show available data arrays or vtk_data[\"varname\"] to retrieve the VTKDataArray for variable varname.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ReadVTK.VTKDataArray","page":"Reference","title":"ReadVTK.VTKDataArray","text":"VTKDataArray{T, N, Format}\n\nHold information about a single VTK data array (cell data or point data).\n\nThe data type is encoded as T, N represents the size of the second dimension in case of multi-dimensonal arrays (or 1 for a vector), and Format encodes in which format the data is stored in the XML file.\n\nThe actual data can be retrieved by calling get_data on the PVTKDataArray object.\n\nSee also: get_data\n\n\n\n\n\n","category":"type"},{"location":"reference/#ReadVTK.VTKDataArray-Tuple{Any, VTKFile}","page":"Reference","title":"ReadVTK.VTKDataArray","text":"VTKDataArray(xml_element, vtk_file)\n\nCreate a lightweight container for a given xml_element and the corresponding vtk_file.\n\nThe xml_element must be of type DataArray and the vtk_file parameter is required for retrieving the actual data. You can pass a VTKDataArray object to get_data to retrieve the actual data.\n\nSee also: get_data\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.VTKFile","page":"Reference","title":"ReadVTK.VTKFile","text":"VTKFile\n\nHold all relevant information about a VTK XML file that has been read in.\n\nFields\n\nfilename: original path to the VTK file that has been read in\nxml_file: object that represents the XML file\nfile_type: currently only \"UnstructuredGrid\" or \"ImageData\" are supported\nversion: VTK XML file format version\nbyte_order: can be LittleEndian or BigEndian and must currently be the same as the system's\ncompressor: can be empty (no compression) or vtkZLibDataCompressor\nappended_data: in case of appended data (see XML documentation), the data is stored here for                  convenient retrieval (otherwise it is empty)\nn_points: number of points in the VTK file\nn_cells: number of cells in the VTK file`\n\n\n\n\n\n","category":"type"},{"location":"reference/#ReadVTK.VTKFile-Tuple{Any}","page":"Reference","title":"ReadVTK.VTKFile","text":"VTKFile(filename)\n\nRead in and parse the VTK XML file specified by its filename.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.VTKPrimitives","page":"Reference","title":"ReadVTK.VTKPrimitives","text":"VTKPrimitives{Connectivity, Offsets}\n\nStore the connectivity and offsets information from the VTK PolyData file as one-dimensional array-like containers. See the VTK file format documentation for information on how to connect the connectivity and offset arrays to the actual geometric points.\n\nYou may use length to determine the number of cells from a VTKPrimitives object.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ReadVTK.get_cell_data-Tuple{PVTKFile}","page":"Reference","title":"ReadVTK.get_cell_data","text":"get_cell_data(pvtk_file::PVTKFile)\n\nRetrieve a lightweight vector with PVTKData objects with the cell data of the given PVTK files. Only numeric data (i.e., DataArray) elements will be read.\n\nSee also: PVTKData, get_cell_data\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_cell_data-Tuple{VTKFile}","page":"Reference","title":"ReadVTK.get_cell_data","text":"get_cell_data(vtk_file::VTKFile)\n\nRetrieve a lightweight VTKData object with the cell data of the given VTK file. Only numeric data (i.e., DataArray) elements will be read.\n\nSee also: VTKData, get_point_data\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_cells-Tuple{VTKFile}","page":"Reference","title":"ReadVTK.get_cells","text":"get_cells(vtk_file)\n\nRetrieve VTK cell information as an object of type VTKCells.\n\nSee also: VTKCells\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_coordinate_data-Tuple{PVTKFile}","page":"Reference","title":"ReadVTK.get_coordinate_data","text":"get_coordinate_data(pvtk_file::PVTKFile)\n\nRetrieve a lightweight {VTKData object with the coordinate data of the given VTK file. Only coordinates of numeric data (i.e., DataArray) elements will be read.\n\nSee also: PVTKData, get_point_data,  get_cell_data\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_coordinate_data-Tuple{VTKFile}","page":"Reference","title":"ReadVTK.get_coordinate_data","text":"get_coordinate_data(vtk_file::VTKFile)\n\nRetrieve a lightweight VTKData object with the coordinate data of the given VTK file. Only coordinates of numeric data (i.e., DataArray) elements will be read.\n\nSee also: VTKData, get_point_data,  get_cell_data\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_coordinates-Tuple{PVTKFile}","page":"Reference","title":"ReadVTK.get_coordinates","text":"get_coordinates(pvtk_file::{VTKFile; x_string=\"x\", y_string=\"y\", z_string=\"z\")\n\nRetrieve VTK coordinate vectors in each direction as a tuple of 1D vectors for a PRectilinearGrid file.\n\nNote that in VTK, points are always stored three-dimensional, even for 1D or 2D files, so you will always retrieve a tuple with three vectors.\n\nSee also: get_cells\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_coordinates-Tuple{VTKFile}","page":"Reference","title":"ReadVTK.get_coordinates","text":"get_coordinates(vtk_file::VTKFile; x_string=\"x\", y_string=\"y\", z_string=\"z\")\n\nRetrieve VTK coordinate vectors in each direction as a tuple of 1D vectors for a RectilinearGrid file.\n\nNote that in VTK, points are always stored three-dimensional, even for 1D or 2D files, so you will always retrieve a tuple with three vectors.\n\nSee also: get_cells\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_data","page":"Reference","title":"ReadVTK.get_data","text":"get_data(data_array::VTKDataArray)\n\nRetrieve actual data from a VTKDataArray as a one- or two-dimensional array-like container.\n\nNote: This function is not type stable but could be - help wanted!\n\n\n\n\n\n","category":"function"},{"location":"reference/#ReadVTK.get_data-Tuple{PVTKDataArray}","page":"Reference","title":"ReadVTK.get_data","text":"get_data(data_array::PVTKDataArray)\n\nRetrieve actual data from a PVTKDataArray as a one- or two-dimensional array-like container.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_data_reshaped-Tuple{PVTKDataArray}","page":"Reference","title":"ReadVTK.get_data_reshaped","text":"Retrieve actual data from a PVTKDataArray and reshapes it as 1D, 2D, or 3D arrays, in case we deal with parallel structured grids. It also puts it in the correct location the the full grid\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_data_reshaped-Union{Tuple{VTKDataArray{T, N, Format} where Format}, Tuple{N}, Tuple{T}} where {T, N}","page":"Reference","title":"ReadVTK.get_data_reshaped","text":"get_data_reshaped(data_array::VTKDataArray; cell_data=false)\n\nRetrieve actual data from a VTKDataArray and reshape it as 1D, 2D, or 3D arrays, in case we deal with structured grids. Note that vectors or tensors will have their components stored in the first dimension of the array. As there is no way to automatically tell from the VTK file format whether it is a tensor, the user has to reshape this manually.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_example_file-Tuple{Any}","page":"Reference","title":"ReadVTK.get_example_file","text":"get_example_file(filename; head=\"main\", output_directory=\".\", force=false)\n\nRetrieve an example file from the ReadVTK_examples repository at commit/branch head and store it in the output_directory. If the file already exists locally, do not download the file again unless force is true. Return the local path to the downloaded file.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_extents","page":"Reference","title":"ReadVTK.get_extents","text":"get_extents(xml_file, min_extent=[0;0;0])\n\nRetrieve the local size of pieces of a structured grid (ImageData, RectilinearGrid). Note that this always returns three dimensions, even if the data is 1D or 2D.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ReadVTK.get_imagedata_dataset-Tuple{PVTKFile}","page":"Reference","title":"ReadVTK.get_imagedata_dataset","text":"get_imagedata_dataset(pvtk_file::PVTKFile)\n\nRetrieve a vector of ImageData datasets from the given PVTKFile file.\n\nSee also: PVTKFile\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_imagedata_dataset-Tuple{VTKFile}","page":"Reference","title":"ReadVTK.get_imagedata_dataset","text":"get_imagedata_dataset(vtk_file::VTKFile)\n\nRetrieve ImageData dataset from the given VTKFile file.\n\nSee also: VTKFile\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_origin-Tuple{PVTKFile}","page":"Reference","title":"ReadVTK.get_origin","text":"get_origin(pvtk_file::PVTKFile)\n\nRetrieve the vector of coordinates of the origin of a uniform grid from the given PVTKFile file.\n\nSee also: PVTKFile\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_origin-Tuple{VTKFile}","page":"Reference","title":"ReadVTK.get_origin","text":"get_origin(vtk_file::VTKFile)\n\nRetrieve the vector of coordinates of the origin of a uniform grid from the given VTKFile file.\n\nSee also: VTKFile\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_point_data-Tuple{PVTKFile}","page":"Reference","title":"ReadVTK.get_point_data","text":"get_point_data(pvtk_file::PVTKFile)\n\nRetrieve a lightweight vector with PVTKData objects with the point data of the given PVTK files. Only numeric data (i.e., DataArray) elements will be read.\n\nSee also: PVTKData, get_cell_data\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_point_data-Tuple{VTKFile}","page":"Reference","title":"ReadVTK.get_point_data","text":"get_point_data(vtk_file::VTKFile)\n\nRetrieve a lightweight VTKData object with the point data of the given VTK file. Only numeric data (i.e., DataArray) elements will be read.\n\nSee also: VTKData, get_cell_data\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_points-Tuple{PVTKFile}","page":"Reference","title":"ReadVTK.get_points","text":"getpoints(vtkfile::PVTKFile)\n\nRetrieve VTK points as a two-dimensional array-like container for a parallel file\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_points-Tuple{VTKFile}","page":"Reference","title":"ReadVTK.get_points","text":"get_points(vtk_file::VTKFile)\n\nRetrieve VTK points as a two-dimensional array-like container.\n\nThe points are stored in dimension × points format. Note that in VTK, points are always stored three-dimensional, even for 1D or 2D files.\n\nSee also: get_cells\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_primitives-Tuple{VTKFile, AbstractString}","page":"Reference","title":"ReadVTK.get_primitives","text":"getprimitives(vtkfile, primitive_type::AbstractString)\n\nRetrieve VTK primitives as an object of type VTKPrimitives. Supported values of primitive type are : \"Verts\", \"Lines\", or \"Polys\".\n\nSee also: VTKPrimitives\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_spacing-Tuple{PVTKFile}","page":"Reference","title":"ReadVTK.get_spacing","text":"get_spacing(pvtk_file::PVTKFile)\n\nRetrieve a vector with the regular increments in each coordinate direction of the uniform grid from the given PVTKFile file.\n\nSee also: PVTKFile\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_spacing-Tuple{VTKFile}","page":"Reference","title":"ReadVTK.get_spacing","text":"get_spacing(vtk_file::VTKFile)\n\nRetrieve a vector with the regular increments in each coordinate direction of the uniform grid from the given VTKFile file.\n\nSee also: VTKFile\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_whole_extent-Tuple{PVTKFile}","page":"Reference","title":"ReadVTK.get_whole_extent","text":"get_whole_extent(pvtk_file::PVTKFile)\n\nRetrieve a vector with the WholeExtent 6-entry vector from the uniform grid PVTKFile file.\n\nSee also: PVTKFile\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_whole_extent-Tuple{VTKFile}","page":"Reference","title":"ReadVTK.get_whole_extent","text":"get_whole_extent(vtk_file::VTKFile)\n\nRetrieve a vector with the WholeExtent 6-entry vector from the uniform grid VTKFile file.\n\nSee also: VTKFile\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.get_wholeextent","page":"Reference","title":"ReadVTK.get_wholeextent","text":"getwholeextent(xmlfile, cell_data=false)\n\nRetrieve the size of a structured grid (ImageData, RectilinearGrid). Note that this always returns three dimensions, even if the data is 1D or 2D.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ReadVTK.isstructured-Tuple{Any}","page":"Reference","title":"ReadVTK.isstructured","text":"isstructured(xml_file)\n\nReturns true if it is a structured grid.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ReadVTK.to_meshcells-Tuple{VTKCells}","page":"Reference","title":"ReadVTK.to_meshcells","text":"to_meshcells(cells::VTKCells)\n\nConvert a VTKCells object, which holds raw point and connectivity data for a number of cells, to a vector of MeshCell objects. The latter can, e.g., be passed to the WriteVTK.jl package for writing new VTK files. See also: VTKCells, VTKBase.MeshCell\n\n\n\n\n\n","category":"method"},{"location":"external/VTKBase/#VTKBase.jl","page":"VTKBase.jl","title":"VTKBase.jl","text":"","category":"section"},{"location":"external/VTKBase/","page":"VTKBase.jl","title":"VTKBase.jl","text":"The VTKBase.jl package contains common definitions used in the WriteVTK.jl and ReadVTK.jl packages.","category":"page"},{"location":"external/VTKBase/#VTK-dataset-types","page":"VTKBase.jl","title":"VTK dataset types","text":"","category":"section"},{"location":"external/VTKBase/","page":"VTKBase.jl","title":"VTKBase.jl","text":"AbstractVTKDataset\nStructuredVTKDataset\nVTKImageData\nVTKRectilinearGrid\nVTKStructuredGrid\nUnstructuredVTKDataset\nVTKPolyData\nVTKUnstructuredGrid","category":"page"},{"location":"external/VTKBase/#VTKBase.AbstractVTKDataset","page":"VTKBase.jl","title":"VTKBase.AbstractVTKDataset","text":"AbstractVTKDataset\n\nAbstract type representing any structured or unstructured VTK dataset.\n\nThe dataset classification is described in the VTK file format specification, page 12.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.StructuredVTKDataset","page":"VTKBase.jl","title":"VTKBase.StructuredVTKDataset","text":"StructuredVTKDataset <: AbstractVTKDataset\n\nAbstract type representing a structured VTK dataset.\n\nSubtypes are VTKImageData, VTKRectilinearGrid and VTKStructuredGrid.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.VTKImageData","page":"VTKBase.jl","title":"VTKBase.VTKImageData","text":"VTKImageData <: StructuredVTKDataset\n\nRepresents the VTK image data format (.vti extension).\n\nThis corresponds to rectangular grids with uniform spacing in all directions.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.VTKRectilinearGrid","page":"VTKBase.jl","title":"VTKBase.VTKRectilinearGrid","text":"VTKRectilinearGrid <: StructuredVTKDataset\n\nRepresents the VTK rectilinear grid format (.vtr extension).\n\nThis corresponds to rectangular grids with non-uniform spacing.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.VTKStructuredGrid","page":"VTKBase.jl","title":"VTKBase.VTKStructuredGrid","text":"VTKStructuredGrid <: StructuredVTKDataset\n\nRepresents the VTK structured grid format (.vts extension).\n\nThis corresponds to curvilinear grids, the most general kind of structured grid.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.UnstructuredVTKDataset","page":"VTKBase.jl","title":"VTKBase.UnstructuredVTKDataset","text":"UnstructuredVTKDataset <: AbstractVTKDataset\n\nAbstract type representing an unstructured VTK dataset.\n\nSubtypes are VTKPolyData and VTKUnstructuredGrid.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.VTKPolyData","page":"VTKBase.jl","title":"VTKBase.VTKPolyData","text":"VTKPolyData <: UnstructuredVTKDataset\n\nRepresents the VTK polydata format (.vtp extension).\n\nThese are unstructured datasets that accept a limited set of cells types, defined in the PolyData module.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.VTKUnstructuredGrid","page":"VTKBase.jl","title":"VTKBase.VTKUnstructuredGrid","text":"VTKUnstructuredGrid <: UnstructuredVTKDataset\n\nRepresents the VTK unstructured format (.vtu extension).\n\nThis is the most general kind of unstructured grid, which accepts all cell types defined in the VTKCellTypes module.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#Field-data-types","page":"VTKBase.jl","title":"Field data types","text":"","category":"section"},{"location":"external/VTKBase/","page":"VTKBase.jl","title":"VTKBase.jl","text":"In VTK, data can either be attached to the geometry (point and cell data), or not (field data).","category":"page"},{"location":"external/VTKBase/","page":"VTKBase.jl","title":"VTKBase.jl","text":"VTKBase.AbstractFieldData\nVTKPointData\nVTKCellData\nVTKFieldData","category":"page"},{"location":"external/VTKBase/#VTKBase.AbstractFieldData","page":"VTKBase.jl","title":"VTKBase.AbstractFieldData","text":"AbstractFieldData\n\nAbstract type representing any kind of dataset.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.VTKPointData","page":"VTKBase.jl","title":"VTKBase.VTKPointData","text":"VTKPointData <: AbstractFieldData\n\nRepresents data that is to be attached to grid points.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.VTKCellData","page":"VTKBase.jl","title":"VTKBase.VTKCellData","text":"VTKCellData <: AbstractFieldData\n\nRepresents data that is to be attached to grid cells.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.VTKFieldData","page":"VTKBase.jl","title":"VTKBase.VTKFieldData","text":"VTKFieldData <: AbstractFieldData\n\nRepresents data that is not attached to the grid geometry.\n\nThis is typically used for lightweight metadata, such as timestep information or strings.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#Cells-in-unstructured-grids","page":"VTKBase.jl","title":"Cells in unstructured grids","text":"","category":"section"},{"location":"external/VTKBase/#General-unstructured-datasets","page":"VTKBase.jl","title":"General unstructured datasets","text":"","category":"section"},{"location":"external/VTKBase/","page":"VTKBase.jl","title":"VTKBase.jl","text":"These are useful when working with general unstructured datasets (.vtu files).","category":"page"},{"location":"external/VTKBase/","page":"VTKBase.jl","title":"VTKBase.jl","text":"VTKCellTypes\nVTKCellTypes.nodes\nVTKBase.AbstractMeshCell\nMeshCell","category":"page"},{"location":"external/VTKBase/#VTKBase.VTKCellTypes","page":"VTKBase.jl","title":"VTKBase.VTKCellTypes","text":"VTKCellTypes\n\nModule defining cell types for unstructured datasets.\n\nDefinitions are adapted from the VTK source code.\n\n\n\n\n\n","category":"module"},{"location":"external/VTKBase/#VTKBase.VTKCellTypes.nodes","page":"VTKBase.jl","title":"VTKBase.VTKCellTypes.nodes","text":"nodes(c::VTKCellTypes)\n\nReturns the number of nodes (or grid points) required by the cell type.\n\nFor instance, this returns 3 for VTK_TRIANGLE.\n\nFor cell types that can take any number of nodes, such as VTK_POLY_LINE, this returns -1.\n\n\n\n\n\n","category":"function"},{"location":"external/VTKBase/#VTKBase.AbstractMeshCell","page":"VTKBase.jl","title":"VTKBase.AbstractMeshCell","text":"AbstractMeshCell\n\nAbstract type specifying a VTK cell.\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#VTKBase.MeshCell","page":"VTKBase.jl","title":"VTKBase.MeshCell","text":"MeshCell <: AbstractMeshCell\n\nSingle cell element in unstructured or polygonal grid.\n\nIt is characterised by a cell type (for instance, VTKCellType.TRIANGLE or PolyData.Strips) and by a connectivity vector determining the points on the grid defining this cell.\n\n\n\nMeshCell(cell_type, connectivity)\n\nDefine a single cell element of an unstructured grid.\n\nThe cell_type argument characterises the type of cell (e.g. vertex, triangle, hexaedron, ...):\n\ncell types for unstructured datasets are defined in the VTKCellTypes\n\nmodule;\n\ncell types for polygonal datasets are defined in the PolyData module.\n\nThe connectivity argument is a vector or tuple containing the indices of the points passed to vtk_grid which define this cell.\n\nExample\n\nDefine a triangular cell passing by points with indices [3, 5, 42].\n\njulia> cell = MeshCell(VTKCellTypes.VTK_TRIANGLE, (3, 5, 42))\nMeshCell{VTKCellType, Tuple{Int64, Int64, Int64}}(VTKCellType(\"VTK_TRIANGLE\", 0x05, 3), (3, 5, 42))\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#Polygonal-datasets","page":"VTKBase.jl","title":"Polygonal datasets","text":"","category":"section"},{"location":"external/VTKBase/","page":"VTKBase.jl","title":"VTKBase.jl","text":"These are useful when working with polygonal datasets (.vtp files).","category":"page"},{"location":"external/VTKBase/","page":"VTKBase.jl","title":"VTKBase.jl","text":"PolyData\nVTKPolyhedron","category":"page"},{"location":"external/VTKBase/#VTKBase.PolyData","page":"VTKBase.jl","title":"VTKBase.PolyData","text":"PolyData\n\nDefines cell types for polygonal datasets.\n\nThe following singleton types are defined:\n\nPolyData.Verts for vertices,\nPolyData.Lines for lines,\nPolyData.Strips for triangular strips,\nPolyData.Polys for polygons.\n\n\n\n\n\n","category":"module"},{"location":"external/VTKBase/#VTKBase.VTKPolyhedron","page":"VTKBase.jl","title":"VTKBase.VTKPolyhedron","text":"VTKPolyhedron <: AbstractMeshCell\n\nRepresents a polyhedron cell in an unstructured grid.\n\nUsing VTKPolyhedron should be preferred to using a MeshCell with a cell type VTKCellTypes.VTK_POLYHEDRON, since the latter cannot hold all the necessary information to describe a polyhedron cell.\n\n\n\nVTKPolyhedron(connectivity, faces...)\n\nConstruct polyhedron cell from connectivity vector (see MeshCell for details) and from a list of polyhedron faces.\n\nExample\n\nCreate a polyhedron with 8 points and 6 faces. This can represent a cube if the 8 points are properly positioned.\n\njulia> cell = VTKPolyhedron(\n           1:8,\n           (1, 4, 3, 2),\n           (1, 5, 8, 4),\n           (5, 6, 7, 8),\n           (6, 2, 3, 7),\n           (1, 2, 6, 5),\n           (3, 4, 8, 7),\n       )\nVTKPolyhedron{UnitRange{Int64}, NTuple{6, NTuple{4, Int64}}}(1:8, ((1, 4, 3, 2), (1, 5, 8, 4), (5, 6, 7, 8), (6, 2, 3, 7), (1, 2, 6, 5), (3, 4, 8, 7)))\n\n\n\n\n\n","category":"type"},{"location":"external/VTKBase/#Constants","page":"VTKBase.jl","title":"Constants","text":"","category":"section"},{"location":"external/VTKBase/","page":"VTKBase.jl","title":"VTKBase.jl","text":"VTKBase.VTKDataType","category":"page"},{"location":"external/VTKBase/#VTKBase.VTKDataType","page":"VTKBase.jl","title":"VTKBase.VTKDataType","text":"VTKDataType\n\nUnion of integer, float and string data types allowed by VTK.\n\n\n\n\n\n","category":"type"},{"location":"#ReadVTK.jl","page":"Home","title":"ReadVTK.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"With ReadVTK.jl you can read in data from VTK XML files in Julia. It aims to complement the excellent package WriteVTK.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note: ReadVTK was mainly motivated by wanting to write proper tests for Trixi2Vtk.jl. A lot of useful features are still missing (see What does not work), and community contributions to improve this package are welcome!","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"First, load the package with","category":"page"},{"location":"","page":"Home","title":"Home","text":"using ReadVTK","category":"page"},{"location":"","page":"Home","title":"Home","text":"Open a VTK file by creating a VTKFile object and passing the filename to the constructor:","category":"page"},{"location":"","page":"Home","title":"Home","text":"vtk = VTKFile(get_example_file(\"celldata_appended_binary_compressed.vtu\"))","category":"page"},{"location":"","page":"Home","title":"Home","text":"To retrieve information about the cell data, use","category":"page"},{"location":"","page":"Home","title":"Home","text":"cell_data = get_cell_data(vtk)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The return object of type VTKCellData allows access to the individual VTKDataArrays using a dictionary-like syntax:","category":"page"},{"location":"","page":"Home","title":"Home","text":"element_ids = cell_data[\"element_ids\"]","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, the actual data can be obtained by executing","category":"page"},{"location":"","page":"Home","title":"Home","text":"data = get_data(element_ids)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Full example including REPL output:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using ReadVTK\n\njulia> vtk = VTKFile(get_example_file(\"celldata_appended_binary_compressed.vtu\"))\nVTKFile(\"celldata_appended_binary_compressed.vtu\", <XMLDocument>, \"UnstructuredGrid\", \"1.0.0\", \"LittleEndian\", \"vtkZLibDataCompressor\", <appended_data>, 4434, 3085)\n\njulia> cell_data = get_cell_data(vtk)\nVTKData()\n\njulia> element_ids = cell_data[\"element_ids\"]\nVTKDataArray(\"element_ids\")\n\njulia> data = get_data(element_ids)\n3085-element reinterpret(Int64, ::Vector{UInt8}):\n    1\n    2\n    3\n    ⋮\n 3083\n 3084\n 3085","category":"page"},{"location":"","page":"Home","title":"Home","text":"After modifications to the read VTK data, one can write back using WriteVTK.jl but must first convert cell objects using to_meshcells.  Continuing from the REPL code above:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using WriteVTK\n\njulia> points = get_points(vtk); cells = to_meshcells(get_cells(vtk));\n\njulia> vtk_grid(\"celldata_appended_binary_compressed_new.vtu\", points, cells) do vtk\n         vtk[\"element_ids\"] = data\n       end\n1-element Vector{String}:\n \"celldata_appended_binary_compressed_new.vtu\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Further example VTK files can be found in the ReadVTK_examples repository.","category":"page"},{"location":"#What-works","page":"Home","title":"What works","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Reading in VTK XML files of type UnstructuredGrid, StructuredGrid, RectilinearGrid,ImageData, PUnstructuredGrid, PStructuredGrid, PRectilinearGrid,PImageData, or PolyData\nExtracting cell or point data\nExtracting point coordinates\nExtracting information about cell types\nOnly for ImageData,PImageData files: get origin, spacing, and extent information\nOnly for RectilinearGrid,PRectiLinearGrid files: get 1D coordinate vectors\nOnly for StructuredGrid,PStructuredGrid files: get coordinate arrays\nReading PolyData files containing vortices, lines, and/or polygons\nReading PVD files\nReading ParaView VTK files that are in-line binary (experimental, only UnstructuredGrid type tested)","category":"page"},{"location":"#What-does-not-work","page":"Home","title":"What does not work","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Reading VTK files not stored in the VTK XML format\nReading VTK files of other type than what is listed under What works above\nMultiblock files\nDifferent byte orders in file and host system\nProbably reading from VTK files that were not created by WriteVTK.jl will fail, specifically since\ncompressed data is assumed to be stored as a single block\nappended data is assumed to be stored as raw\nheader_type is hardcoded to UInt64\nExtracting primitives from PolyData files other than vortices, lines, and/or polygons\nLikely anything else that is not specifically mentioned under What works","category":"page"},{"location":"#Development","page":"Home","title":"Development","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Helpful resources for working with (i.e., reading and writing) VTK XML files:","category":"page"},{"location":"","page":"Home","title":"Home","text":"VTK file format documentation (incomplete!) as a PDF\nVTK XML formats wiki article\nBlog post on encoding binary data\nMailing list message on encoding binary data","category":"page"},{"location":"","page":"Home","title":"Home","text":"We use JuliaFormatter.jl to keep a consistent code formatting. If you have installed JuliaFormatter.jl, just run","category":"page"},{"location":"","page":"Home","title":"Home","text":"using JuliaFormatter; format(\".\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"in the top-level directory of ReadVTK.jl to update the formatting.","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ReadVTK is maintained by the Trixi authors. Its principal developers are Michael Schlottke-Lakemper (University of Stuttgart, Germany) and Hendrik Ranocha (Johannes Gutenberg University Mainz, Germany).","category":"page"},{"location":"","page":"Home","title":"Home","text":"Further contributions to ReadVTK have been made by the following people:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Jorge Pérez Zerpa (Universidad de la República, Uruguay)\nOndřej Kincl (Charles University, Czech Republic)\nBoris Kaus (Johannes-Gutenberg University Mainz, Germany)\nMatthew Whisenant (University of Tennessee, Knoxville)","category":"page"},{"location":"#License-and-contributing","page":"Home","title":"License and contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ReadVTK is licensed under the MIT license (see License). Since ReadVTK is an open-source project, we are very happy to accept contributions from the community. Please refer to Contributing for more details. To get in touch with the developers, join us on Trixi's Slack workspace or create an issue.","category":"page"},{"location":"#Acknowledgments","page":"Home","title":"Acknowledgments","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package would not exist without the nice work of Juan Ignacio Polanco and his cleanly written and well-documented package WriteVTK.jl.","category":"page"}]
}
