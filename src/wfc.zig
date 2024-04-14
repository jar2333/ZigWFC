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

pub const SquareTile = struct {
    xpos: u32,
    xneg: u32,
    ypos: u32,
    yneg: u32
};

const SquareDimensions = struct {
    x: u32,
    y: u32
};

pub const CubeTile = struct {
    xpos: u32,
    xneg: u32,
    ypos: u32,
    yneg: u32,
    zpos: u32,
    zneg: u32
};

const CubeDimensions = struct {
    x: u32,
    y: u32,
    z: u32
};

pub const WFCError = error{
    InvalidGridSize
};

pub fn Solver(comptime TileT: type) type {
    // Comptime constants
    const DimensionT: type = switch (TileT) {
        SquareTile => SquareDimensions,
        CubeTile => CubeDimensions,
        else => @compileError("Tile type " ++ @typeName(TileT) ++ " not supported.")
    };
    const dim: comptime_int = @typeInfo(DimensionT).Struct.fields.len;

    return struct {
        const Self = @This();

        allocator: std.mem.Allocator = undefined,
        tileset: []const TileT = undefined,

        pub fn init(alloc: std.mem.Allocator, tiles: []const TileT) Self {
            // Create an internal data structure for representing adjacencies using Tile data
            
            // Return solver
            return .{
                .allocator = alloc,
                .tileset = tiles,
            };
        }

        pub fn solve(_: *Self, grid: []usize, dimensions: DimensionT) !void {
            // Check if provided dimensions fit into provided grid
            var prod: u32 = 1;
            inline for (std.meta.fields(DimensionT)) |f| {
                prod *= @as(u32, @field(dimensions, f.name));
            }
            if (prod != grid.len) {
                return WFCError.InvalidGridSize;
            }

            
            return;
        }

        const TileIndex = usize;
        const GridIndex = usize;

        const SideIndex = u8;

        const BitsetT = std.bit_set.StaticBitSet(5);

        //
        // ==============
        // = Variables
        // ==============
        // tiles: []TileT
        // Array containing each tile's adjacency data. Index with usize.
        //
        // grid: []TileIndex, grid.len == prod(dimensions)  
        // Flat array that contains indexes to Tiles array. Index with usize. 
        // This index represents a grid position. Write a helper that will help yield the x,y,z position coordinates.
        // 
        // possibilities: []BitsetT, possibilities.len == grid.len
        // Array which corresponds to grid array (grid[i] -> neighbors[i]). It is an array of bitsets.
        // The bitset at possibilities[i] corresponds to the set of possible tiles at tile index i.
        // To be modified at each iteration of WFC in one of two ways:
        // 1) Collapse: Make only one tile be possible at this index. Set exactly one bit to 1, the rest to 0.
        // 2) Propagate: Propagate the results of a collapse at the ith grid index using DFS. This means until stack is empty:
        //   2a) Find all neighbor grid indeces and their associated position in neighbors[i] (make an array n, n.len==neighbors[i].len)
        //   2b) Using the bitset at neighbors[i][j] and the one at possibilities[i] to update the one at possibilities[n[j]]
        //   2c) Add the neighbors key n[j] to the stack if the propagation caused a change to possibilities[n[j]]
        // 
        // neighbors: [][num_sides]const BitsetT, neighbors.len == tiles.len, bitset.capacity == tiles.len 
        // Array which corresponds to tiles array (tiles[i] -> neighbors[i]). It is an array of arrays of bitsets. 
        // The bitset[i][j] corresponds to the set of possible tiles relative to tile i at neighbor tile j. The number of possible tiles is the length of tiles array.
        // Basically, each tile has num_sides neighbors, and each of those has a possibility space (hence the bitset) relative to the current tile. 
        // Created by processing the initial tiles array, and by taking account the grid shape. The possibilities are fixed hence const.
        // To be queried to propagate the results of collapsing a tile.
        //
        fn processInitialConstraints() void {}
        fn isCollapsed() bool {return false;}
        fn isContradiction() bool {return false;}
        fn iterate() void {}
        fn getMinEntropyCoordinates() usize {return 0;}
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

    var solver = Solver(SquareTile).init(allocator, &tiles);

    var grid = try allocator.alloc(usize, 100);
    defer allocator.free(grid);

    // Test that it errors out when provided dimensions do not match grid 
    const err = try solver.solve(grid, .{.x = 5, .y = 10});
    testing.expect(err == WFCError.InvalidGridSize);

    try solver.solve(grid, .{.x = 10, .y = 10});
}