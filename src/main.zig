const std = @import("std");
const testing = std.testing;

const wfc = @import("./wfc.zig");

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    const allocator = arena.allocator();

    const tiles = [_]wfc.SquareTile{
        wfc.SquareTile{.xpos = 1, .xneg = 1, .ypos = 1, .yneg = 1}
    };

    var prng = std.rand.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.os.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    const rand = prng.random();
    var solver = try wfc.Solver(wfc.SquareTile).init(allocator, &tiles, rand);
    defer solver.deinit();

    var grid = try allocator.alloc(usize, 100);
    defer allocator.free(grid);

    // Test that it errors out when provided dimensions do not match grid 
    // try solver.solve(grid, .{.x = 5, .y = 10});
    // testing.expectError(WFCError.InvalidGridSize, err);

    try solver.solve(grid, .{.x = 10, .y = 10});

}