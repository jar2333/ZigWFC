const std = @import("std");

// --------------------------
// - The proposed interface:
// --------------------------
//
// Internal grid types: 
// - Square grid
// - Cube grid
//
// Input:
// - 2D array / slice of slices -> square grid
// - 3d array / slice of slices of slices -> cube grid
//
const Solver = struct {

    const allocator: *std.mem.Allocator = undefined;
    const Self = @This();

    pub fn init(alloc: *std.mem.Allocator) Self {
        return .{
            .allocator = alloc,
        };
    }

    pub fn solve(_: Self) void {
        return;
    }

};