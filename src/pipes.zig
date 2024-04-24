const std = @import("std");
const unicode = @import("std").unicode;
const io = @import("std").io;

const wfc = @import("./wfc.zig");

const CLIError = error {
    NoValueProvided,
    InvalidArgumentName,
    InvalidArgumentNameLength,
    UnsupportedArgumentOption,
    InsufficientArguments,
};

const WFCConfig = struct {
    seed: ?u64,
    tiles: []const wfc.SquareTile,
    width: usize,
    height: usize,
};

const sample_tiles = [_]wfc.SquareTile{
    // ∘
    wfc.SquareTile{.xpos = 1, .ypos = 1, .xneg = 1, .yneg = 1},
    // ─
    wfc.SquareTile{.xpos = 2, .ypos = 1, .xneg = 2, .yneg = 1},
    // │
    wfc.SquareTile{.xpos = 1, .ypos = 2, .xneg = 1, .yneg = 2},
    // ┼
    wfc.SquareTile{.xpos = 2, .ypos = 2, .xneg = 2, .yneg = 2},
    // ┌
    wfc.SquareTile{.xpos = 2, .ypos = 2, .xneg = 1, .yneg = 1},
    // ┘
    wfc.SquareTile{.xpos = 1, .ypos = 1, .xneg = 2, .yneg = 2},
    // └
    wfc.SquareTile{.xpos = 2, .ypos = 1, .xneg = 1, .yneg = 2},
    // ┐
    wfc.SquareTile{.xpos = 1, .ypos = 2, .xneg = 2, .yneg = 1},
    // ┬ 
    wfc.SquareTile{.xpos = 2, .ypos = 2, .xneg = 2, .yneg = 1},
    // ┤
    wfc.SquareTile{.xpos = 1, .ypos = 2, .xneg = 2, .yneg = 2},
    // ├
    wfc.SquareTile{.xpos = 2, .ypos = 2, .xneg = 1, .yneg = 2},
    // ┴
    wfc.SquareTile{.xpos = 2, .ypos = 1, .xneg = 2, .yneg = 2},
};

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
        else => "◦"
    });
}

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    const allocator = arena.allocator();

    // ==================================
    // = Process command line arguments
    // ==================================
    const cfg_opt = try processCommandLineArgs(allocator);
    
    if (cfg_opt) |cfg| {
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
}

fn processCommandLineArgs(allocator: std.mem.Allocator) !?WFCConfig {
    var argsIterator = try std.process.ArgIterator.initWithAllocator(allocator);
    defer argsIterator.deinit();

    // Skip executable name
    _ = argsIterator.next();

    var cfg = WFCConfig{
        .seed = null, 
        .tiles = &sample_tiles,
        .width = 25,
        .height = 10,
    };

    while (argsIterator.next()) |arg| {
        // Check for argument formatting
        if (arg[0] != '-') {
            return CLIError.InvalidArgumentName;
        }
        if (arg.len > 2) {
            return CLIError.InvalidArgumentNameLength;
        }

        // Process arguments with no value
        if (arg[1] == 'h') {
            try printHelp();
            return null;
        }

        // Process arguments with 1 value
        const value: [:0]const u8 = argsIterator.next() orelse return CLIError.NoValueProvided;

        switch (arg[1]) {
            'x' => cfg.width  = try std.fmt.parseInt(usize, value, 10),
            'y' => cfg.height = try std.fmt.parseInt(usize, value, 10),
            's' => cfg.seed   = try std.fmt.parseInt(u64, value, 10),
            else => return CLIError.UnsupportedArgumentOption
        }
    }

    return cfg;
}

const helpMessage = 
    \\
    \\This demo prints a solved grid with the sample tiles, using the provided width, height, and RNG seed.
    \\
    \\Tiles:
    \\ ◦ ─ │ ┌ ┘ └ ┐ ┬ ┤ ├ ┴ ┼
    \\
    \\Args:
    \\  -x: Specifies the grid width. Default value of 25.
    \\  -y: Specifies the grid height. Default value of 10.
    \\  -s: Specifies the seed for the RNG. Default is a random value.
    \\  -h: Prints this message.
    \\
;

fn printHelp() !void {
    const w = io.getStdOut().writer();
    try w.print(helpMessage, .{});
}

fn runWFC(allocator: std.mem.Allocator, cfg: WFCConfig, grid: []usize) !void {
    const seed = cfg.seed;
    const tiles = cfg.tiles;
    const width = cfg.width;
    const height = cfg.height;

    var prng = std.rand.DefaultPrng.init(seed orelse blk:{
        var s: u64 = undefined;
        try std.posix.getrandom(std.mem.asBytes(&s));
        break :blk s;
    });
    const rand = prng.random();

    var solver = try wfc.Solver(wfc.SquareTile).init(allocator, tiles, rand);
    defer solver.deinit();

    try solver.solve(grid, .{.x = width, .y = height});
}

fn printSquareGrid(grid: []const usize, width: usize, height: usize) !void {
    std.debug.assert(grid.len == width*height);
    const w = io.getStdOut().writer();

    try w.print("{} x {}:\n", .{width, height});

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