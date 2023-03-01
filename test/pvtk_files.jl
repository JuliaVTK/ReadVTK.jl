# test parallel VTK files
using WriteVTK, Test, ReadVTK

# (1) Generate parallel structured input files

# Global grid
xs_global = range(0, 2; length = 15)
ys_global = range(-1, 1; length = 12)
zs_global = range(0, 1; length = 4)
v_global =  range(0,1,length=3)

extents = [
    ( 1:10,  1:5, 1:4),  # process 1
    (10:15,  1:5, 1:4),  # process 2
    ( 1:10, 5:12, 1:4),  # process 3
    (10:15, 5:12, 1:4),  # process 4
]

saved_files = Vector{Vector{String}}(undef, 4)  # files saved by each "process"

# Write *.pvti file
for part = 1:4
    is, js, ks = extents[part]  # local indices
    xs, ys, zs = xs_global[is], ys_global[js], zs_global[ks]  # local grid
    xs_c, ys_c, zs_c = xs[1:end-1], ys[1:end-1], zs[1:end-1]

    saved_files[part] = pvtk_grid(
            "fields", xs, ys, zs;
            part = part, extents = extents,
        ) do pvtk
        pvtk["Temperature"] = [x + 2y + 3z for x ∈ xs, y ∈ ys, z ∈ zs]
        pvtk["Velocity"] = [x + 2y + 3z + v for v ∈ v_global, x ∈ xs, y ∈ ys, z ∈ zs]
        pvtk["Pressure"] = [x + 2y + 3z for x ∈ xs_c, y ∈ ys_c, z ∈ zs_c]
        pvtk["Phase"] = [trunc(Int64, x*15) for x ∈ xs_c, y ∈ ys_c, z ∈ zs_c]
    end
end

# write *.pvtr file
for part = 1:4
  is, js, ks = extents[part]  # local indices
  xs, ys, zs = xs_global[is], ys_global[js], zs_global[ks]  # local grid
  xs_c, ys_c, zs_c = xs[1:end-1], ys[1:end-1], zs[1:end-1]
  saved_files[part] = pvtk_grid(
          "fields", Vector(xs), Vector(ys), Vector(zs);
          part = part, extents = extents,
      ) do pvtk
      pvtk["Temperature"] = [x + 2y + 3z for x ∈ xs, y ∈ ys, z ∈ zs]
      pvtk["Velocity"] = [x + 2y + 3z + v for v ∈ v_global, x ∈ xs, y ∈ ys, z ∈ zs]
      pvtk["Pressure"] = [x + 2y + 3z for x ∈ xs_c, y ∈ ys_c, z ∈ zs_c]
      pvtk["Phase"] = [trunc(Int64, x*15) for x ∈ xs_c, y ∈ ys_c, z ∈ zs_c]
  end
end

T_global = [x + 2y + 3z for x ∈ xs_global, y ∈ ys_global, z ∈ zs_global]
V_global = [x + 2y + 3z + v for v ∈ v_global, x ∈ xs_global, y ∈ ys_global, z ∈ zs_global]
P_global = [x + 2y + 3z for x ∈ xs_global[1:end-1], y ∈ ys_global[1:end-1], z ∈ zs_global[1:end-1]]
Phase_global = [trunc(Int64, x*15) for x ∈ xs_global[1:end-1], y ∈ ys_global[1:end-1], z ∈ zs_global[1:end-1]]
        

# Unstructured grid
all_data = [
    # Process 1
    (
        points = rand(3, 5),  # 5 points on process 1
        cells = [             # 2 cells  on process 1
            MeshCell(VTKCellTypes.VTK_TRIANGLE, [1, 4, 2]),
            MeshCell(VTKCellTypes.VTK_QUAD, [2, 4, 3, 5]),
        ],
    ),

    # Process 2
    (
        points = rand(3, 4),  # 4 points on process 2
        cells = [             # 1 cell   on process 2
            MeshCell(VTKCellTypes.VTK_QUAD, [1, 2, 3, 4]),
        ]
    ),
]

saved_files = Vector{Vector{String}}(undef, 2)  # files saved by each "process"

for part = 1:2
    data = all_data[part]
    saved_files[part] = pvtk_grid(
            "simulation", data.points, data.cells;
            part = part, nparts = 2,
        ) do pvtk
        pvtk["Pressure"] = sum(data.points; dims = 1)
    end
end




# (2) Read back files

# a) RectilinearGrid file
pvtk = PVTKFile("fields.pvtr")

coords_read = get_coordinates(pvtk)
@test Vector(xs_global) == coords_read[1]
@test Vector(ys_global) == coords_read[2]
@test Vector(zs_global) == coords_read[3]

# Extract data
point_data = get_point_data(pvtk)
T_read = get_data_reshaped(point_data["Temperature"])
V_read = get_data_reshaped(point_data["Velocity"])
@test  T_global == T_read
@test  V_global == V_read

cell_data = get_cell_data(pvtk)
P_read = get_data_reshaped(cell_data["Pressure"], cell_data=true)
Phase_read = get_data_reshaped(cell_data["Phase"], cell_data=true)
@test P_global == P_read
@test Phase_global == Phase_read

# ImageData file
pvtk = PVTKFile("fields.pvti")
origin = get_origin(pvtk)
@test origin[1] == xs_global[1]
@test origin[2] == ys_global[1]
@test origin[3] == zs_global[1]

# Extract data
point_data = get_point_data(pvtk)
T_read = get_data_reshaped(point_data["Temperature"])
V_read = get_data_reshaped(point_data["Velocity"])
@test  T_global == T_read
@test  V_global == V_read

cell_data = get_cell_data(pvtk)
P_read = get_data_reshaped(cell_data["Pressure"], cell_data=true)
Phase_read = get_data_reshaped(cell_data["Phase"], cell_data=true)
@test  P_global == P_read
@test Phase_global == Phase_read


# c) PVTU files
pvtk = PVTKFile("simulation.pvtu")
points = get_points(pvtk)
@test all_data[1].points == points[1]
@test all_data[2].points == points[2]

point_data = get_point_data(pvtk)
P_read = get_data(point_data["Pressure"])
@test sum(all_data[1].points,dims=1)[:] == P_read[1]
@test sum(all_data[2].points,dims=1)[:] == P_read[2]


