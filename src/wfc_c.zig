const std = @import("std");

const wfc = @import("wfc.zig");

const SquareSolver = wfc.Solver(wfc.SquareTile, .{});
const CubeSolver = wfc.Solver(wfc.CubeTile, .{});
const WFCError = wfc.WFCError;

// ====================
// = C API
// ====================

// --------------
// - Error codes
// --------------

export const WFC_OK: c_int = 0;
export const WFC_CONTRADICTION: c_int = 1;
export const WFC_INVALID_GRID_SIZE: c_int = 2;
export const WFC_TOO_MANY_TILES: c_int = 3;
export const WFC_OUT_OF_MEMORY: c_int = 4;

// ---------------
// - Square grids
// ---------------
const wfc_SquareTile = wfc.SquareTile;

// Returns a void pointer. It is null if an error ocurred, an error code is set in this case.
export fn wfc_initSquareGridSolver(tiles: [*]wfc_SquareTile, tiles_size: usize, seed: u64, err_code: *c_int) ?*anyopaque {
    const allocator = std.heap.c_allocator;

    const solver = allocator.create(SquareSolver) catch {
        err_code.* = WFC_OUT_OF_MEMORY;
        return null;
    };

    var prng = std.rand.DefaultPrng.init(seed);
    const rand = prng.random();
    solver.* = SquareSolver.init(allocator, tiles[0..tiles_size], rand) catch |err| {
        switch (err) {
            WFCError.TooManyTiles => err_code.* = WFC_TOO_MANY_TILES,
            else => err_code.* = WFC_OUT_OF_MEMORY, 
        }
        return null;
    };

    err_code.* = WFC_OK;
    return solver;
}

export fn wfc_freeSquareGridSolver(solver: *anyopaque) void {
    const solver_ptr: *SquareSolver = @ptrCast(@alignCast(solver));
    solver_ptr.deinit();
}

export fn wfc_solveSquareGrid(solver: *anyopaque, grid: [*]u8, grid_size: usize, width: usize, height: usize) c_int {
    const solver_ptr: *SquareSolver = @ptrCast(@alignCast(solver));

    solver_ptr.solve(grid[0..grid_size], .{.x = width, .y = height}) catch |err| {
        switch (err) {
            WFCError.Contradiction => return WFC_CONTRADICTION,
            WFCError.InvalidGridSize => return WFC_INVALID_GRID_SIZE,
            else => return WFC_OUT_OF_MEMORY, 
        }
    };

    return WFC_OK;
}
