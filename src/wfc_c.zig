const std = @import("std");

const wfc = @import("wfc.zig");

const SquareSolver = wfc.Solver(wfc_SquareTile, .{});
const WFCError = wfc.WFCError;

// ====================
// = C API
// ====================

const wfc_SquareTile = wfc.SquareTile;

// Returns a void pointer. It is null if an error ocurred, an error code is set in this case.
export fn wfc_initSquareGridSolver(tiles: [*]wfc_SquareTile, tiles_size: usize, seed: u64, err_code: *c_int) ?*anyopaque {
    const allocator = std.heap.c_allocator;

    const solver = allocator.create(SquareSolver) catch {
        err_code.* = -10;
        return null;
    };

    var prng = std.rand.DefaultPrng.init(seed);
    const rand = prng.random();
    solver.* = SquareSolver.init(allocator, tiles[0..tiles_size], rand) catch |err| {
        switch (err) {
            WFCError.TooManyTiles => err_code.* = -3,
            else =>                  err_code.* = -10, 
        }
        return null;
    };

    return solver;
}

export fn wfc_freeSquareGridSolver(solver: *anyopaque) void {
    const solver_ptr: *SquareSolver = @ptrCast(@alignCast(solver));
    solver_ptr.deinit();

    // Is this implicit free idiomatic?
    // const allocator = std.heap.c_allocator;
    // allocator.free(solver_ptr);
}

export fn wfc_solveSquareGrid(solver: *anyopaque, grid: [*]u8, grid_size: usize, width: usize, height: usize) c_int {
    const solver_ptr: *SquareSolver = @ptrCast(@alignCast(solver));

    solver_ptr.solve(grid[0..grid_size], .{.x = width, .y = height}) catch |err| {
        switch (err) {
            WFCError.Contradiction => return -1,
            WFCError.InvalidGridSize => return -2,
            else => return -10, 
        }
    };

    return 0;
}
