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
    const DimensionT: type = switch (TileT) {
        SquareTile => SquareDimensions,
        CubeTile => CubeDimensions,
        else => @compileError("Tile type " ++ @typeName(TileT) ++ " not supported.")
    };

    return struct {
        const Self = @This();

        allocator: std.mem.Allocator = undefined,
        tileset: []const TileT = undefined,

        pub fn init(alloc: std.mem.Allocator, tiles: []const TileT) Self {
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