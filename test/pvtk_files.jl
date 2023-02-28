# test parallel VTK files
using WriteVTK

# (1) Generate parallel structured input files

# Global grid
xs_global = range(0, 2; length = 15)
ys_global = range(-1, 1; length = 12)
zs_global = range(0, 1; length = 4)

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
    saved_files[part] = pvtk_grid(
            "fields", xs, ys, zs;
            part = part, extents = extents,
        ) do pvtk
        pvtk["Temperature"] = [x + 2y + 3z for x ∈ xs, y ∈ ys, z ∈ zs]
    end
end

# write *.pvtr file
for part = 1:4
  is, js, ks = extents[part]  # local indices
  xs, ys, zs = xs_global[is], ys_global[js], zs_global[ks]  # local grid
  saved_files[part] = pvtk_grid(
          "fields", Vector(xs), Vector(ys), Vector(zs);
          part = part, extents = extents,
      ) do pvtk
      pvtk["Temperature"] = [x + 2y + 3z for x ∈ xs, y ∈ ys, z ∈ zs]
  end
end


# (2) Read back files




