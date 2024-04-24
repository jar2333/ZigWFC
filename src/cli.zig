const std = @import("std");
const unicode = @import("std").unicode;
const io = @import("std").io;

const wfc = @import("./wfc.zig");

const CLIError = error {
    NoValueProvided,
    InvalidArgumentName,
    InsufficientArguments,
    Error
};

const WFCConfig = struct {
    seed: ?u64,
    tiles: []const wfc.SquareTile,
    width: usize,
    height: usize,
};

const sample_tiles = [_]wfc.SquareTile{
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

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    const allocator = arena.allocator();

    // ==================================
    // = Process command line arguments
    // ==================================
    const cfg = processCommandLineArgs(allocator);

    // =============
    // = Run WFC
    // =============
    const grid = try allocator.alloc(usize, cfg.width*cfg.height);
    defer allocator.free(grid);

    try runWFC(allocator, cfg, grid);

    // ===============
    // = Print grid
    // ===============
    try printSquareGrid(grid, cfg.width, cfg.height);
}

fn processCommandLineArgs(allocator: std.mem.Allocator) !WFCConfig {
    var argsIterator = try std.process.ArgIterator.initWithAllocator(allocator);
    defer argsIterator.deinit();

    // Skip executable name
    _ = argsIterator.next();

    var cfg = WFCConfig{
        .seed = null, 
        .tiles = sample_tiles,
        .width = 10,
        .height = 10,
    };

    while (argsIterator.next()) |arg| {
        if (arg[0] != '-') {
            return CLIError.InvalidArgumentName;
        }
        if (arg.len == 0 or arg.len > 1) {
            return CLIError.InvalidArgumentName;
        }

        const value: [:0]const u8 = argsIterator.next() orelse return CLIError.NoValueProvided;

        switch (arg[1]) {
            'w' => cfg.width  = try std.fmt.parseInt(usize, value, 10),
            'h' => cfg.height = try std.fmt.parseInt(usize, value, 10),
            's' => cfg.seed   = try std.fmt.parseInt(u64, value, 10),
            't' => cfg.tiles  = try parseTiles(allocator, value),
            else => return CLIError.InvalidArgumentName
        }
    }

    return cfg;
}

// Parse tile data from a file
fn parseTiles(allocator: std.mem.Allocator, _: [:0]const u8) ![]wfc.SquareTiles {


    return try allocator.alloc(wfc.SquareTile, 0);
}

fn runWFC(allocator: std.mem.Allocator, cfg: WFCConfig, grid: []usize) !void {
    const seed_opt = cfg.seed;
    const tiles = cfg.tiles;
    const width = cfg.width;
    const height = cfg.height;

    var prng = if (seed_opt) |seed| {
        std.rand.DefaultPrng.init(seed);
    }
    else {
        var seed: u64 = undefined;
        try std.posix.getrandom(std.mem.asBytes(&seed));
        std.rand.DefaultPrng.init(seed);
    };

    const rand = prng.random();

    var solver = try wfc.Solver(wfc.SquareTile).init(allocator, tiles, rand);
    defer solver.deinit();

    try solver.solve(grid, .{.x = width, .y = height});
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