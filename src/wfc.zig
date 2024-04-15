const std = @import("std");

// --------------------------
// - The proposed interface:
// --------------------------
//
// Internal grid types: 
// - Square grid
// - Cube grid
// - Hexagonal grid
//
// Input:
// - A const array of structs specifying adjacency data for tiles.
// - A mutable (flat, multidimensional) array/slice that will be populated with the indexes of the tiles

// Implementation detail:
// The Tile types are packed structs. This allows treating them like 
// circular buffers through @bitCast, to efficiently find opposite sides.
// i.e. The opposite side to tile[i] is tile[i+(num_sides/2) % num_sides]

// The canonical side ordering for the packed structs / ring buffers is:
// 1) xpos (right)
// 2) ypos (up)
// 3) zpos (?)
// 4) xneg (left)
// 5) yneg (down)
// 6) zneg (?)

const LabelT = u32;

pub const SquareTile = packed struct {
    xpos: LabelT = 0,
    ypos: LabelT = 0,
    xneg: LabelT = 0,
    yneg: LabelT = 0
};

const SquareDimensions = packed struct {
    x: u32,
    y: u32
};

pub const CubeTile = packed struct {
    xpos: LabelT = 0,
    ypos: LabelT = 0,
    zpos: LabelT = 0,
    xneg: LabelT = 0,
    yneg: LabelT = 0,
    zneg: LabelT = 0
};

const CubeDimensions = packed struct {
    x: u32,
    y: u32,
    z: u32
};

pub const WFCError = error{
    InvalidGridSize,
    Contradiction
};

// O(tiles.len*num_sides) complexity adjacency extraction algorithm:
// Preliminary breakdown: 
// Given a tile and a side of that tile, we want to query all other tiles that connect on the opposite side.
// 
// The input to the solver is an array of structs, each of which contain labels that denote tile adjacency. 
// Two sides of different tiles can be adjacent if and only if their labels are equal.
//
// Implementation detail:
// We can efficiently query a side's label and find the opposite side if we can represent the tile labels struct as 
// a contiguous array of memory, to be treated as a ring buffer.
//
// To this end, we use Zig packed structs in the API. This allows the user to comfortably set the labels for a tile's sides.
// However, it also allows us to guarantee the correct memory layout needed to later use the struct as a ring buffer.
// We can @bitCast the tile struct to an array of labels. A side can be indexed by k, and the label queried as buffer[k].
// The index of the opposite side is thus (k+(n/2))%n, where n is the number of sides in a tile.

// To solve the initial problem, we can use a bucketing algorithm, following this reasoning:
// For k in 0..n:
//   Given a tile with label A at side k, it connects to tiles which have label A at side (k+(n/2))%n.

// Thus, we can iterate over every tile, and add them to n buckets. Each bucket corresponds to a particular (label, side_index) pair. 
// For each of a tile's side indeces k with label buffer[k], we add the tile's index i to bucket[buffer[k]][(k+(n/2))%n]. 
// 
// After iterating over all tiles once, each possible (label, side_index) pair should index a bucket with every tile whose opposite side index contains that label 
//
// In pseudocode:
// 
// bucket = hashmap<labeltype, list<tile_index>[n]>
//
// for i in  0..tiles.len:
//   t = tiles[i]
//   buffer = bitcast<labeltype[n]>(t)
//   for k in 0..n:
//     label = buffer[k]
//     opposite = (k+(n/2))%n
//     bucket[label][opposite].add(i)

pub fn Solver(comptime TileT: type) type {
    // Comptime constants
    const DimensionT: type = switch (TileT) {
        SquareTile => SquareDimensions,
        CubeTile => CubeDimensions,
        else => @compileError("Tile type " ++ @typeName(TileT) ++ " not supported.")
    };
    const dim: comptime_int = @typeInfo(DimensionT).Struct.fields.len;
    const num_sides: comptime_int = @typeInfo(TileT).Struct.fields.len;

    return struct {
        const Self = @This();

        const TileIndex = usize;
        const GridIndex = usize;

        const SideIndex = usize;

        const BitsetT = std.bit_set.DynamicBitSet;

        allocator: std.mem.Allocator = undefined,

        // ==============
        // = Variables
        // ==============
        //
        // tileset: []TileT
        // Array containing each tile's adjacency data. Index with usize (TileIndex).
        //
        // grid: []TileIndex, grid.len == prod(dimensions)  
        // Flat array that contains indexes to Tiles array. Index with usize (GridIndex). 
        // This index represents a grid position. Write a helper that will help yield the x,y,z position coordinates.
        // 
        // possibilities: []BitsetT, possibilities.len == grid.len
        // Array which corresponds to grid array (grid[i] -> neighbors[i]). It is an array of bitsets.
        // The bitset at possibilities[i] corresponds to the set of possible tiles at tile index i.
        // To be modified at each iteration of WFC in one of two ways:
        // 1) Collapse: Make only one tile be possible at this index. Set exactly one bit to 1, the rest to 0.
        // 2) Propagate: Propagate the results of a collapse at the ith grid index using DFS. This means until stack is empty:
        //   2a) Find all neighbor grid indeces ordered by their associated bitset in neighbors[i] (make an array n: [num_sides]GridIndex, n.len==neighbors[i].len)
        //   2b) Using the bitset at neighbors[i][j] and the one at possibilities[i] to update the one at possibilities[n[j]]
        //   2c) Add the neighbors key n[j] to the stack if the propagation caused a change to possibilities[n[j]]
        // 
        // neighbors: [][num_sides]const BitsetT, neighbors.len == tiles.len, bitset.capacity == tiles.len 
        // Array which corresponds to tiles array (tiles[i] -> neighbors[i]). It is an array of arrays of bitsets. 
        // The bitset[i][j] corresponds to the set of possible tiles relative to tile i at neighbor tile j. The number of possible tiles is the length of tiles array.
        // Basically, each tile has num_sides neighbors, and each of those has a possibility space (hence the bitset) relative to the current tile. 
        // Created by processing the initial tiles array, and by taking account the grid shape. The possibilities are fixed hence const.
        // To be queried to propagate the results of collapsing a tile.
        
        // At init
        tileset: []const TileT = undefined,
        neighbors: []const[num_sides]BitsetT = undefined,

        // At solve
        dimensions: DimensionT = undefined,
        grid: []TileIndex = undefined,
        possibilities: []BitsetT = undefined,
        
        pub fn init(alloc: std.mem.Allocator, tiles: []const TileT) !Self {
            // Initialize neighbors array and bitsets in neighbors array
            const neighbors: []const[num_sides]BitsetT = try alloc.alloc([num_sides]BitsetT, tiles.len);
            for (0..tiles.len) |i| {
                for (0..num_sides) |k| {
                    neighbors[i][k] = try BitsetT.initEmpty(alloc, tiles.len);
                }
            }

            // Use bucketing algorithm to map all (label, side_index) pairs to a list of possible adjacent tiles
            // Room for optimization: Instead of list, use an array bounded by tiles.len and a capacity. If tiles is known at comptime, no dynamic allocation needed.
            const Bucket = std.ArrayList(TileIndex);
            const HashMap = std.array_hash_map.AutoArrayHashMap(LabelT, [num_sides]Bucket);
            var buckets = HashMap.init(alloc);
            defer buckets.deinit();

            for (0..tiles.len) |i| {
              const buffer: [num_sides]LabelT = @bitCast(tiles[i]);
              for (0..num_sides) |k| {
                const label = buffer[k];
                const opposite_k = (k+(num_sides/2))%num_sides;

                // Lazy initialization of all num_sides adjacency lists
                if (!buckets.contains(label)) {
                    try buckets.put(label, [num_sides]Bucket{
                        Bucket.initCapacity(alloc, tiles.len),
                        Bucket.initCapacity(alloc, tiles.len),
                        Bucket.initCapacity(alloc, tiles.len),
                        Bucket.initCapacity(alloc, tiles.len),
                    });
                    for (buckets.getPtr(label).?.*.items) |b| {
                        defer b.deinit();
                    }
                }
                // Append the tile index i to the adjacency list
                const ptr: ?*[num_sides]Bucket = buckets.getPtr(label);
                ptr.?.*[opposite_k].append(i);
              }
            }

            // Populate neighbors array by setting bitsets
            for (tiles, 0..) |tile, i| {
                // neighbors: []const[num_sides]BitsetT => neighbors[tile_index][side_index].set(neighbor_index)
                // buckets: std.array_hash_map.AutoArrayHashMap(LabelT, [num_sides]std.ArrayList(TileIndex)) => buckets.getPtr(label).?.*[opposite_side_index][]
                const buffer: [num_sides]LabelT = @bitCast(tile);
                for (0..num_sides) |k| {
                    const label = buffer[k];

                    // For each tile index in the adjacency list, set the corresponding bit in the bitset
                    const ptr: ?*[num_sides]Bucket = buckets.getPtr(label);
                    for (ptr.?.*[k].items) |j| {
                        neighbors[i][k].set(j);
                    }
                }
            }

            // Return solver
            return .{
                .allocator = alloc,
                .tileset = tiles,
                .neighbor = neighbors,
            };
        }

        pub fn deinit(_: *Self) void {}

        pub fn solve(self: *Self, grid: []usize, dimensions: DimensionT) !void {
            // Check if provided dimensions fit into provided grid
            var size: u32 = 1;
            inline for (std.meta.fields(DimensionT)) |f| {
                const v = @as(u32, @field(dimensions, f.name));
                std.debug.assert(v != 0);
                size *= v;
            }
            if (size != grid.len) {
                return WFCError.InvalidGridSize;
            }

            // Solving algorithm:
            self.grid = grid;
            self.dimensions = dimensions;

            while (!isCollapsed()) {
                iterate();
            }

            //can be moved down the call stack if required by algorithm logic, later
            if (isContradiction()) {
                return WFCError.Contradiction;
            }
        }

        fn isCollapsed(_: *Self) bool {
            return false;
        }
        
        fn isContradiction(_: *Self) bool {
            return false;
        }

        fn iterate(_: *Self) void {

        }

        fn getMinEntropyCoordinates(_: *Self) GridIndex {
            return 0;
        }

        fn collapseAt(_: *Self, _: GridIndex) void {

        }

        fn propagate(_: *Self, _: GridIndex) void {

        }

        //returns true if neighbor's possible tiles decrease
        fn propagateAt(_: *Self, current: GridIndex, neighbor: GridIndex) bool {
            return false;
        }

        
        // Row major indexing formulas:
        // 2D
        // i = x + width*y =>
        // x = i % width;
        // y = i / width;
        // 3D: 
        // i = x + width*y + width*height*z =>
        // x = i % width;
        // y = (i / width)%height;
        // z = i / (width*height)
        // Note: Do i use a *[dim]GridIndex or a []GridIndex with an assert(len == dim) ? 
        // TODO: Add bounds checking!!! Hell, I updated the neighbors type to be a slice of optionals
        fn getNeighbors(self: *Self, neighbors: []?GridIndex, p: GridIndex) void {
            std.debug.assert(neighbors.len == dim);

            // Gets easily indexable width, height(, depth)
            const buffer: [dim]u32 = @bitCast(self.dimensions);

            // Treat on case by case basis due to special logic for each grid type
            switch (TileT) {
                SquareTile => {
                    const width = buffer[0];

                    const x = p % width;
                    const y = p / width;

                    neighbors[0] = (x+1) + width*y; // xpos neighbor
                    neighbors[1] = x + width*(y+1); // ypos neighbor
                    neighbors[2] = (x-1) + width*y; // xneg neighbor
                    neighbors[3] = x + width*(y-1); // yneg neighbor
                },
                CubeTile => {
                    const width  = buffer[0];
                    const height = buffer[1];
                    
                    const x = p % width;
                    const y = (p / width)%height;
                    const z = p / (width*height);

                    neighbors[0] = (x+1) + width*y + width*height*z; // xpos neighbor
                    neighbors[1] = x + width*(y+1) + width*height*z; // ypos neighbor
                    neighbors[3] = x + width*y + width*height*(z+1); // zpos neighbor
                    neighbors[4] = (x-1) + width*y + width*height*z; // xpos neighbor
                    neighbors[5] = x + width*(y-1) + width*height*z; // ypos neighbor
                    neighbors[6] = x + width*y + width*height*(z-1); // zpos neighbor
                },
                else => @compileError("Tile type " ++ @typeName(TileT) ++ " not yet supported.")
            }

        }

        fn collapseRandom(_: *Self, tiles: *const BitsetT) TileIndex {
            tiles.set(0);
            return 0;
        }

        fn getAdjacencies(self: *Self, p: GridIndex, k: SideIndex) *BitsetT {
            return &self.neighbors[p][k];
        }

        fn getPossibleTiles(self: *Self, p: GridIndex) *BitsetT {
            return &self.possibilities[p];
        }

    };
}

const testing = std.testing;
const test_allocator = std.testing.allocator;

test "basic solver test" {
    const allocator = test_allocator;

    const tiles = [_]SquareTile{
        SquareTile{.xpos = 0, .xneg = 0, .ypos = 0, .yneg = 0}
    };

    var solver = try Solver(SquareTile).init(allocator, &tiles);

    var grid = try allocator.alloc(usize, 100);
    defer allocator.free(grid);

    // Test that it errors out when provided dimensions do not match grid 
    const err = try solver.solve(grid, .{.x = 5, .y = 10});
    testing.expect(err == WFCError.InvalidGridSize);

    try solver.solve(grid, .{.x = 10, .y = 10});
}