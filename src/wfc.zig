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
// - An array of structs specifying adjacency data for tiles.
// - A (flat) array that will be populated with the indexes of the tiles
// - An enum value that determines the type of grid/tiles used (maybe can be omitted through reflection on 1st arg) 

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
    if (TileT != SquareTile and TileT != CubeTile) {
        @compileError("Tile type " ++ @typeName(TileT) ++ " not supported.");
    }
    return struct {
        const Self = @This();
        const DimensionT: type = switch (TileT) {
            SquareTile => SquareDimensions,
            CubeTile => CubeDimensions,
            else => @compileError("Tile type " ++ @typeName(TileT) ++ " not supported.")
        };

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

test "basic add functionality" {
    const allocator = test_allocator;

    const tiles = [_]SquareTile{
        SquareTile{.xpos = 0, .xneg = 0, .ypos = 0, .yneg = 0}
    };

    var solver = Solver(SquareTile).init(allocator, &tiles);

    var grid = try allocator.alloc(usize, 100);
    defer allocator.free(grid);

    try solver.solve(grid, .{.x = 10, .y = 10});
}