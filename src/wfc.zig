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

pub const SquareTile = packed struct {
    xpos: u32 = 0,
    ypos: u32 = 0,
    xneg: u32 = 0,
    yneg: u32 = 0
};

const SquareDimensions = struct {
    x: u32,
    y: u32
};

pub const CubeTile = packed struct {
    xpos: u32 = 0,
    ypos: u32 = 0,
    zpos: u32 = 0,
    xneg: u32 = 0,
    yneg: u32 = 0,
    zneg: u32 = 0
};

const CubeDimensions = struct {
    x: u32,
    y: u32,
    z: u32
};

pub const WFCError = error{
    InvalidGridSize,
    Contradiction
};

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
        grid: []TileIndex = undefined,
        possibilities: []BitsetT = undefined,
        
        pub fn init(alloc: std.mem.Allocator, tiles: []const TileT) !Self {
            // Initialize neighbors array
            const neighbors: []const[num_sides]BitsetT = try alloc.alloc([num_sides]BitsetT, tiles.len);

            // Init bitsets in neighbors array
            for (0..tiles.len) |i| {
                for (0..num_sides) |j| {
                    neighbors[i][j] = BitsetT.initEmpty(alloc, tiles.len);
                }
            }

            // Set the corresponding bits by iterating over all pairs of tiles, and noting when one side of one corresponds to the opposite side in the other
            // This is a naive O(n^2) algorithm with redundancy, try finding a more efficient one. (bitset operations, or something else)
            for (tiles, 0..) |tile, i| {
                for (tiles, 0..) |other, j| {
                    for (0..num_sides) |k| {
                        // Reinterpret tile packed struct as array of u32 to map side to opposite side mathematically
                        const tile_arr: [num_sides]u32 = @bitCast(tile);
                        const other_arr: [num_sides]u32 = @bitCast(other);

                        // Modify neighbors[i][k] and neighbors[j][ok] if a match is found]
                        const ok = @mod(k+(num_sides/2), num_sides);
                        if (tile_arr[k] == other_arr[ok]) {
                            neighbors[i][k].set(j);
                            neighbors[j][ok].set(i);
                        }
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

        pub fn solve(self: *Self, grid: []usize, dimensions: DimensionT) !void {
            // Check if provided dimensions fit into provided grid
            var size: u32 = 1;
            inline for (std.meta.fields(DimensionT)) |f| {
                size *= @as(u32, @field(dimensions, f.name));
            }
            if (size != grid.len) {
                return WFCError.InvalidGridSize;
            }

            // Solving algorithm:
            self.grid = grid;

            while (!isCollapsed()) {
                iterate();
            }

            //can be moved down the call stack if required by algorithm logic, later
            if (isContradiction()) {
                return WFCError.Contradiction;
            }
        }

        fn processInitialConstraints() void {}
        fn isCollapsed() bool {return false;}
        fn isContradiction() bool {return false;}
        fn iterate() void {}
        fn getMinEntropyCoordinates() GridIndex {return 0;}
        fn collapseAt(p: GridIndex) void {}
        fn propagate(p: GridIndex) void {}
        fn propagateAt(current: GridIndex, neighbor: GridIndex) bool {return false;} //returns true if neighbor's possible tiles decrease
        fn getNeighbors(neighbors: [dim]GridIndex, p: GridIndex) void {}

        fn collapseRandom(tiles: *const BitsetT) TileIndex {return 0;}

        fn getAdjacencies(k: GridIndex, d: SideIndex) BitsetT {return BitsetT.init();}

        fn getInitialPossibleTiles(p: GridIndex) BitsetT {}
        fn getPossibleTiles(p: GridIndex) *BitsetT {}

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