const std = @import("std");
const unicode = @import("std").unicode;
const testing = std.testing;

const io = @import("std").io;


const wfc = @import("./wfc.zig");

pub fn main() !void {
    // const allocator = arena.allocator();
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer switch (gpa.deinit()) {
        .ok => std.debug.print("no leaked memory :)\n", .{}),
        .leak => std.debug.print("leaked memory :(\n", .{})
    };

    const allocator = gpa.allocator();

    const tiles = [_]wfc.SquareTile{
        // o
        wfc.SquareTile{.xpos = 1, .ypos = 1, .xneg = 1, .yneg = 1},

        // -
        wfc.SquareTile{.xpos = 2, .ypos = 1, .xneg = 2, .yneg = 1},
        
        // |
        wfc.SquareTile{.xpos = 1, .ypos = 2, .xneg = 1, .yneg = 2},

        // +
        wfc.SquareTile{.xpos = 2, .ypos = 2, .xneg = 2, .yneg = 2},

        // ┌
        wfc.SquareTile{.xpos = 2, .ypos = 2, .xneg = 1, .yneg = 1},

        // ┘
        wfc.SquareTile{.xpos = 1, .ypos = 1, .xneg = 2, .yneg = 2},

        // └
        wfc.SquareTile{.xpos = 2, .ypos = 1, .xneg = 1, .yneg = 2},

        // ┐
        wfc.SquareTile{.xpos = 1, .ypos = 2, .xneg = 2, .yneg = 1},

        wfc.SquareTile{.xpos = 2, .ypos = 2, .xneg = 2, .yneg = 1},

        wfc.SquareTile{.xpos = 1, .ypos = 2, .xneg = 2, .yneg = 2},

        wfc.SquareTile{.xpos = 2, .ypos = 2, .xneg = 1, .yneg = 2},

        wfc.SquareTile{.xpos = 2, .ypos = 1, .xneg = 2, .yneg = 2},
    };

    var prng = std.rand.DefaultPrng.init(blk: {
        var seed: u64 = undefined;
        try std.posix.getrandom(std.mem.asBytes(&seed));
        break :blk seed;
    });
    const rand = prng.random();
    var solver = try wfc.Solver(wfc.SquareTile).init(allocator, &tiles, rand);
    defer solver.deinit();

    const grid = try allocator.alloc(usize, 250);
    for (grid) |*p| {
        p.* = 0;
    }

    defer allocator.free(grid);
    // Test that it errors out when provided dimensions do not match grid 
    // try solver.solve(grid, .{.x = 5, .y = 10});
    // testing.expectError(WFCError.InvalidGridSize, err);

    try solver.solve(grid, .{.x = 25, .y = 10});

    try printSquareGrid(grid, 25, 10);

}

fn printSquareGrid(grid: []const usize, width: usize, height: usize) !void {
    std.debug.assert(grid.len == width*height);
    const w = io.getStdOut().writer();

    try w.print("Printing {} x {} grid:\n", .{width, height});

    var x: usize = 0;
    var y: usize = 0;
    while (y < height): (y += 1) {
        while (x < width): (x += 1) {
            try w.print("{u}", .{try translateUnicode(grid[x+y*width])});
        }
        x = 0;
        try w.writeByte('\n');
    }

    try w.writeByte('\n');
}

fn translateUnicode(index: usize) !u21 {
    return try unicode.utf8Decode(switch (index) {
        1 => "─",
        2 => "│",
        3 => "┼",
        4 => "┌",
        5 => "┘",
        6 => "└",
        7 => "┐",
        8 => "┬",
        9 => "┤",
        10 => "├",
        11 => "┴",
        else => "∘"
    });
}