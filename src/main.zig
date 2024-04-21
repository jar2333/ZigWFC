const std = @import("std");
const testing = std.testing;

const wfc = @import("./wfc.zig");

pub fn main() !void {
    // var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    // defer arena.deinit();

    // const allocator = arena.allocator();
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer switch (gpa.deinit()) {
        .ok => std.debug.print("no leaked memory :)\n", .{}),
        .leak => std.debug.print("leaked memory :(\n", .{})
    };

    const allocator = gpa.allocator();

    const tiles = [_]wfc.SquareTile{
        wfc.SquareTile{.xpos = 1, .ypos = 1, .xneg = 1, .yneg = 1},
        wfc.SquareTile{.xpos = 1, .ypos = 2, .xneg = 3, .yneg = 2},
        wfc.SquareTile{.xpos = 3, .ypos = 3, .xneg = 3, .yneg = 3},
    };

    var prng = std.rand.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.posix.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    const rand = prng.random();
    var solver = try wfc.Solver(wfc.SquareTile).init(allocator, &tiles, rand);
    defer solver.deinit();

    const grid = try allocator.alloc(usize, 100);
    for (grid) |*p| {
        p.* = 0;
    }

    defer allocator.free(grid);
    // Test that it errors out when provided dimensions do not match grid 
    // try solver.solve(grid, .{.x = 5, .y = 10});
    // testing.expectError(WFCError.InvalidGridSize, err);

    try solver.solve(grid, .{.x = 10, .y = 10});

    printSquareGrid(grid, 10, 10);

}

fn printSquareGrid(grid: []const usize, width: usize, height: usize) void {
    std.debug.assert(grid.len == width*height);

    var x: usize = 0;
    var y: usize = 0;
    while (y < height): (y += 1) {
        while (x < width): (x += 1) {
            std.debug.print("{}", .{grid[x+y*width]});
            if (x != width-1) {
                std.debug.print("-", .{});
            }
        }
        x = 0;
        std.debug.print("\n", .{});
    }
}