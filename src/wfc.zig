const std = @import("std");

// --------------------------
// - The proposed interface:
// --------------------------
//         z--------+
//        /|       /|
//       / |      / |
//      o--|-----x  |
//      |  +-----|--+
//      | /      | /
//      |/       |/
//      y--------+
//
// Internal grid types: 
// - Square grid
// - Cube grid
// - Hexagonal grid
//
// Input:
// - A const array of structs specifying adjacency data for tiles.
// - A mutable (flat, multidimensional) array/slice that will be populated with the indexes of the tiles

// Implementation detail:
// The Tile types are packed structs. This allows treating them like 
// circular buffers through @bitCast, to efficiently find opposite sides.
// i.e. The opposite side to tile[i] is tile[i+(num_sides/2) % num_sides]

// The canonical side ordering for the packed structs / ring buffers is:
// 1) xpos (right)
// 2) ypos (up)
// 3) zpos (?)
// 4) xneg (left)
// 5) yneg (down)
// 6) zneg (?)

const LabelT = u8;

pub const SquareTile = packed struct {
    xpos: LabelT = 0,
    ypos: LabelT = 0,
    xneg: LabelT = 0,
    yneg: LabelT = 0
};

const SquareDimensions = packed struct {
    x: usize,
    y: usize
};

pub const CubeTile = packed struct {
    xpos: LabelT = 0,
    ypos: LabelT = 0,
    zpos: LabelT = 0,
    xneg: LabelT = 0,
    yneg: LabelT = 0,
    zneg: LabelT = 0
};

const CubeDimensions = packed struct {
    x: usize,
    y: usize,
    z: usize
};

pub const WFCError = error{
    InvalidGridSize,
    Contradiction,
    TooManyTiles,
};

// O(tiles.len*num_sides) complexity adjacency extraction algorithm:
// Preliminary breakdown: 
// Given a tile and a side of that tile, we want to query all other tiles that connect on the opposite side.
// 
// The input to the solver is an array of structs, each of which contain labels that denote tile adjacency. 
// Two sides of different tiles can be adjacent if and only if their labels are equal.
//
// Implementation detail:
// We can efficiently query a side's label and find the opposite side if we can represent the tile labels struct as 
// a contiguous array of memory, to be treated as a ring buffer.
//
// To this end, we use Zig packed structs in the API. This allows the user to comfortably set the labels for a tile's sides.
// However, it also allows us to guarantee the correct memory layout needed to later use the struct as a ring buffer.
// We can @bitCast the tile struct to an array of labels. A side can be indexed by k, and the label queried as buffer[k].
// The index of the opposite side is thus (k+(n/2))%n, where n is the number of sides in a tile.

// To solve the initial problem, we can use a bucketing algorithm, following this reasoning:
// For k in 0..n:
//   Given a tile with label A at side k, it connects to tiles which have label A at side (k+(n/2))%n.

// Thus, we can iterate over every tile, and add them to n buckets. Each bucket corresponds to a particular (label, side_index) pair. 
// For each of a tile's side indeces k with label buffer[k], we add the tile's index i to bucket[buffer[k]][(k+(n/2))%n]. 
// 
// After iterating over all tiles once, each possible (label, side_index) pair should index a bucket with every tile whose opposite side index contains that label 
//
// In pseudocode:
// 
// bucket = hashmap<labeltype, list<tile_index>[n]>
//
// for i in  0..tiles.len:
//   t = tiles[i]
//   buffer = bitcast<labeltype[n]>(t)
//   for k in 0..n:
//     label = buffer[k]
//     opposite = (k+(n/2))%n
//     bucket[label][opposite].add(i)

// Options allow specifying comptime values to minimize dynamic allocation, or enabling opt-in features.
const SolverOptions = struct {
    enable_callbacks: bool = false,
};

pub fn Solver(comptime TileT: type, comptime _: SolverOptions) type {
    // Comptime constants
    const DimensionT: type = switch (TileT) {
        SquareTile => SquareDimensions,
        CubeTile => CubeDimensions,
        else => @compileError("Tile type " ++ @typeName(TileT) ++ " not supported.")
    };
    const dim: comptime_int = @typeInfo(DimensionT).Struct.fields.len;
    const num_sides: comptime_int = @typeInfo(TileT).Struct.fields.len;

    return struct {
        const Self = @This();

        const TileIndex = u8;
        const GridIndex = usize;

        const MAX_TILES_LENGTH = 1 << @bitSizeOf(TileIndex);

        const BitsetT = std.bit_set.DynamicBitSet;

        // ===================
        // = Member Variables
        // ===================
        //
        // tileset: []TileT
        // Array containing each tile's adjacency data. Index with usize (TileIndex).
        //
        // grid: []TileIndex, grid.len == prod(dimensions)  
        // Flat array that contains indexes to Tiles array. Index with usize (GridIndex). 
        // For SquareTile and CubeTile tilesets, the grid is interpreted as a row-major square/cube grid with given dimensions.
        // This index represents a grid position. Write a helper that will help yield the x,y,z position coordinates.
        // 
        // possibilities: []BitsetT, possibilities.len == grid.len
        // Array which corresponds to grid array (grid[i] -> adjacencies[i]). It is an array of bitsets.
        // The bitset at possibilities[i] corresponds to the set of possible tiles at tile index i.
        // To be modified at each iteration of WFC in one of two ways:
        // 1) Collapse: Make only one tile be possible at this index. Set exactly one bit to 1, the rest to 0.
        // 2) Propagate: Propagate the results of a collapse at the ith grid index using DFS. This means until stack is empty:
        //   2a) Find all neighbor grid indeces ordered by their associated bitset in adjacencies[i] (make an array n: [num_sides]GridIndex, n.len==adjacencies[i].len)
        //   2b) Using the bitset at adjacencies[i][j] and the one at possibilities[i] to update the one at possibilities[n[j]]
        //   2c) Add the adjacencies key n[j] to the stack if the propagation caused a change to possibilities[n[j]]
        // 
        // adjacencies: [][num_sides]const BitsetT, adjacencies.len == tiles.len, bitset.capacity == tiles.len 
        // Array which corresponds to tiles array (tiles[i] -> adjacencies[i]). It is an array of arrays of bitsets. 
        // The bitset[i][j] corresponds to the set of possible tiles relative to tile i at neighbor tile j. The number of possible tiles is the length of tiles array.
        // Basically, each tile has num_sides adjacencies, and each of those has a possibility space (hence the bitset) relative to the current tile. 
        // Created by processing the initial tiles array, and by taking account the grid shape. The possibilities are fixed hence const.
        // To be queried to propagate the results of collapsing a tile.
        
        // At init
        allocator: std.mem.Allocator = undefined,
        tileset: []const TileT = undefined,
        adjacencies: [][num_sides]BitsetT = undefined,
        rand: std.rand.Random = undefined,

        // At solve
        possibilities: []BitsetT = undefined,
        dimensions: DimensionT = undefined,
        grid: []TileIndex = undefined,
        collapsed_count: usize = 0,
        

        /// ==================
        /// = Initialization
        /// ==================

        pub fn init(alloc: std.mem.Allocator, tiles: []const TileT, rand: std.rand.Random) !Self {
            if (tiles.len > MAX_TILES_LENGTH) {
                return WFCError.TooManyTiles;
            }

            // Initialize adjacencies array and bitsets in adjacencies array
            var adjacencies: [][num_sides]BitsetT = try alloc.alloc([num_sides]BitsetT, tiles.len);
            for (0..tiles.len) |i| {
                for (0..num_sides) |k| {
                    adjacencies[i][k] = try BitsetT.initEmpty(alloc, tiles.len);
                }
            }

            try populateAdjacencies(alloc, tiles, adjacencies);
            // std.debug.print("Populated adjacencies:\n", .{});
            // printAdjacencies(adjacencies);

            // Return solver
            return .{
                .allocator = alloc,
                .tileset = tiles,
                .adjacencies = adjacencies,
                .rand = rand
            };
        }

        // Time complexity: O(tiles.len*num_sides)
        fn populateAdjacencies(alloc: std.mem.Allocator, tiles: []const TileT, adjacencies: [][num_sides]BitsetT) !void {
            // Use bucketing algorithm to map all (side_index, label) pairs to a list of possible adjacent tiles
            // NOTE: Instead of list, we use an array bounded by tiles.len and a capacity. If tiles is known at comptime, no dynamic allocation needed.
            const Bucket = std.ArrayListUnmanaged(TileIndex);
            const HashMap = std.AutoArrayHashMap(LabelT, Bucket);

            var buckets: [num_sides]HashMap = [_]HashMap{HashMap.init(alloc)} ** num_sides;
            defer for (&buckets) |*hash| {
                hash.deinit();
            };

            // For every tile, iterate through its sides, and insert the tile to the bucket corresponding to the opposite side and same label
            for (0..tiles.len) |tile_index| {
              const buffer: [num_sides]LabelT = @bitCast(tiles[tile_index]);

              for (0..num_sides) |k| {
                const label: LabelT = buffer[k];
                const opposite_k = (k+(num_sides/2))%num_sides;
                std.debug.assert(0 <= opposite_k and opposite_k <= num_sides-1);

                const hash: *HashMap = &buckets[opposite_k];
                var v = try hash.getOrPut(label);

                // Lazy intitialization: Initialize the bucket at the given side opposite_k if not already there
                if (!v.found_existing) {
                    v.value_ptr.* = try Bucket.initCapacity(alloc, tiles.len);
                }

                // Append the tile index i to the adjacency list
                const new_item_ptr = v.value_ptr.addOneAssumeCapacity();
                new_item_ptr.* = @truncate(tile_index);
              }
            }

            // Deinitialize every bucket as cleanup
            defer for (&buckets) |*hash| {
                for (hash.values()) |*v| {
                    v.deinit(alloc);
                }
            };

            // Populate adjacencies array by setting bitsets
            for (tiles, 0..) |tile, tile_index| {
                // adjacencies: []const[num_sides]BitsetT => adjacencies[tile_index][side_index].set(adjacent_tile_index)
                // buckets: [num_sides]std.AutoHashMap(LabelT, Bucket) => buckets[side_index].getPtr(label).?
                const buffer: [num_sides]LabelT = @bitCast(tile);

                for (0..num_sides) |k| {
                    const label = buffer[k];

                    // For each tile index in the adjacency list, set the corresponding bit in the bitset
                    const opt: ?*Bucket = buckets[k].getPtr(label);
                    if (opt) |ptr| {
                        for (ptr.items) |adjacent_tile_index| {
                            adjacencies[tile_index][k].set(adjacent_tile_index);
                        }
                    }
                }
            }
        }

        /// ==================
        /// = Public API
        /// ==================

        pub fn deinit(self: *Self) void {
            for (0..self.adjacencies.len) |i| {
                for (0..num_sides) |k| {
                    self.adjacencies[i][k].deinit();
                }
            }
            self.allocator.free(self.adjacencies);
        }

        pub fn solve(self: *Self, grid: []u8, dimensions: DimensionT) !void {
            // Check if provided dimensions fit into provided grid
            var size: usize = 1;
            inline for (std.meta.fields(DimensionT)) |f| {
                const v = @as(usize, @field(dimensions, f.name));
                if (v == 0) {
                    return WFCError.InvalidGridSize;
                }
                size *= v;
            }
            if (size != grid.len) {
                return WFCError.InvalidGridSize;
            }
            std.debug.assert(size == grid.len);

            // Solving algorithm setup:
            self.grid = grid;
            self.dimensions = dimensions;
            self.collapsed_count = 0;

            // Initialize possibilities
            self.possibilities = try self.allocator.alloc(BitsetT, grid.len);
            defer self.allocator.free(self.possibilities);

            for (self.possibilities) |*b| {
                b.* = try BitsetT.initFull(self.allocator, self.tileset.len);
            }
            defer for (self.possibilities) |*b| {
                b.deinit();
            };

            // ----------------
            // - Wfc algorithm:
            // ----------------
            while (!self.isCollapsed()) {
                const p: GridIndex = self.getMinEntropyCoordinates();
                try self.collapseAt(p);
                try self.propagate(p);
            }

            // Set all the positions to their solved values
            // Should not panic, since we confirmed all positions have at least 1 possibility
            for (grid, 0..) |*p, i| {
                p.* = @truncate(self.possibilities[i].findFirstSet().?);
            }

        }

        /// ========================
        /// = Private Helpers: WFC
        /// ========================

        // Check if all positions have at most 1 possibility
        fn isCollapsed(self: *Self) bool {
            return self.collapsed_count == self.grid.len;
        }
        
        // Gets the coordinates for the first grid position with the least number of possibilities
        // NOTE: Implementation detail: among all positions with same entropy, the one with lowest index is chosen.
        // NOTE: Naive linear search, can use a memoized result from propagation instead 
        fn getMinEntropyCoordinates(self: *Self) GridIndex {
            var min_entropy_position: GridIndex = 0;
            var min_entropy: usize = self.tileset.len;
            
            for (self.possibilities, 0..) |*p, i|{
                if (p.count() < min_entropy and p.count() > 1) {
                    min_entropy_position = i;
                    min_entropy = p.count();
                }
            }

            return min_entropy_position;
        }

        // NOTE: If alternate collapse behaviors are later supported, modify this function
        fn collapseAt(self: *Self, p: GridIndex) !void {
            // std.debug.print("collapsing at {}\n", .{p});
            const b: *BitsetT = &self.possibilities[p];
            try self.collapseRandom(b);
            self.collapsed_count += 1;
        }

        // DFS traversal propagating the effects of a collapse to rest of grid
        fn propagate(self: *Self, p: GridIndex) !void {
            var stack = std.ArrayList(GridIndex).init(self.allocator);
            defer stack.deinit();

            try stack.append(p);
            while (stack.items.len > 0) {
                const cur = stack.pop();

                var neighbors: [num_sides]?GridIndex = [_]?GridIndex{null} ** num_sides;
                self.getNeighbors(neighbors[0..num_sides], cur);

                for (neighbors, 0..) |opt, k| {
                    if (opt) |n| {
                        if (try self.propagateAt(cur, n, k)) {
                            try stack.append(n);
                        }
                    }
                }
            }
        }

        // Utilize bit operations on the two possibilities bitsets and the adjacency bitset to decrease neighbor count monotonically
        // Explanation:
        //
        // We need to find the set of all allowed tiles for the neighbor grid position at the kth side of the current grid position.
        // To do so:
        // For each tile index i allowed in the possibility space indexed by current (i.e. each index i set to 1 in the bitset getPossibleTiles(current)), 
        // we find the allowed tile indeces for the neighbor at kth side of the ith tile (i.e. the bitset getAdjacencies(i, k))
        // Then, we take the union of all those allowed tiles.
        //
        // After that, we take the intersection of this set with the set of already allowed tiles at the neighbor.
        // Since it is an intersection, the count of allowed tiles decreases monotonically.
        // TODO: Rephrase explanation
        // NOTE: Returns true if neighbor's number of possible tiles decreased, false otherwise
        // NOTE: Following convention of rest of module, k is the index to the side of current tile that is adjacent to neighbor tile
        fn propagateAt(self: *Self, current: GridIndex, neighbor: GridIndex, k: usize) !bool {
            // std.debug.print("Propagating to neighbor {} at side {} from current {}\n", .{neighbor, k, current});
            var neighbor_tiles: *BitsetT = &self.possibilities[neighbor];
            var current_tiles: *BitsetT  = &self.possibilities[current];

            const initial_amount: usize = neighbor_tiles.count();

            var allowed: BitsetT = try BitsetT.initEmpty(self.allocator, self.tileset.len);
            defer allowed.deinit();

            // Efficient iterator using @clz intrinsic
            var iterator = current_tiles.iterator(.{});
            while (iterator.next()) |i| {
                const tile_neighbors: *const BitsetT = self.getAdjacencies(i, k);
                allowed.setUnion(tile_neighbors.*); // efficient bitset union
            }

            neighbor_tiles.setIntersection(allowed); // efficient bitset intersection

            // TODO: Fiddle with the most efficient way to do this part, it is the hot path
            if (neighbor_tiles.count() == initial_amount) {
                return false;
            }

            switch (neighbor_tiles.count()) {
                0 => return WFCError.Contradiction,
                1 => {
                    self.collapsed_count += 1; 
                    return true;
                },
                else => return true
            }
        }
        
        // Row major indexing formulas:
        // 2D
        // i = x + width*y =>
        // x = i % width;
        // y = i / width;
        // 3D: 
        // i = x + width*y + width*height*z =>
        // x = i % width;
        // y = (i / width)%height;
        // z = i / (width*height)
        fn getNeighbors(self: *Self, neighbors: *[num_sides]?GridIndex, p: GridIndex) void {
            const n = self.grid.len;

            // Gets easily indexable width, height(, depth)
            const buffer: [dim]usize = @bitCast(self.dimensions);

            // Treat on case by case basis due to special logic for each grid type
            switch (TileT) {
                SquareTile => {
                    const width = buffer[0];
                    // p = x + width*y =>
                    // x = p % width
                    // y = p / width

                    // Unsigned bounds checking
                    if (p + 1 <= n-1) {
                        neighbors[0] = p + 1;     // xpos neighbor = (x+1) + width*y
                    }
                    if (p + width <= n-1) {
                        neighbors[1] = p + width; // ypos neighbor = x + width*(y+1)
                    }
                    if (p >= 1) { // p - 1 >= 0
                        neighbors[2] = p - 1;     // xneg neighbor = (x-1) + width*y
                    }
                    if (p >= width) { // p - width >= 0
                        neighbors[3] = p - width; // yneg neighbor = x + width*(y-1)
                    }
                },
                CubeTile => {
                    const width  = buffer[0];
                    const height = buffer[1];
                    // p = x + width*y + width*height*z =>
                    // x = p % width
                    // y = (p / width)%height
                    // z = p / (width*height)

                    // Unsigned bounds checking
                    if (p + 1 <= n-1) {
                        neighbors[0] = p + 1;            // xpos neighbor = (x+1) + width*y + width*height*z
                    }
                    if (p + width <= n-1) {
                        neighbors[1] = p + width;        // ypos neighbor = x + width*(y+1) + width*height*z
                    }
                    if (p + width*height <= n-1) {
                        neighbors[2] = p + width*height; // zpos neighbor = x + width*y + width*height*(z+1)
                    }
                    if (p >= 1) { // p - 1 >= 0
                        neighbors[3] = p - 1;            // xneg neighbor = (x-1) + width*y + width*height*z
                    }
                    if (p >= width) { // p - width >= 0
                        neighbors[4] = p - width;        // yneg neighbor = x + width*(y-1) + width*height*z
                    }
                    if (p >= width*height) { // p - width*height >= 0
                        neighbors[5] = p - width*height; // zneg neighbor = x + width*y + width*height*(z-1)
                    }
                },
                else => @compileError("Tile type " ++ @typeName(TileT) ++ " not yet supported.")
            }
        }

        // We need to get indeces of set positions in the bitset, and uniform randomly choose among them.
        // NOTE: Since the absolute max number of tiles is relatively small when TileIndex == u8, we use a fixed buffer allocator :)
        // TODO: Make more efficient using tiles.unsetAll() for static bitsets
        fn collapseRandom(self: *Self, tiles: *BitsetT) !void {
            comptime std.debug.assert(TileIndex == u8);
            var buffer: [MAX_TILES_LENGTH*@sizeOf(u8)]u8 = undefined;
            var fba = std.heap.FixedBufferAllocator.init(&buffer);
            const alloc = fba.allocator();

            var indeces = try std.ArrayListUnmanaged(u8).initCapacity(alloc, tiles.capacity());
            defer indeces.deinit(alloc);

            // Get all set bit indeces and store them, then unset
            var iterator = tiles.iterator(.{});
            while (iterator.next()) |j| {
                const new_item_ptr = indeces.addOneAssumeCapacity();
                new_item_ptr.* = @truncate(j);

                tiles.unset(j);                
            }

            // Set the tile which corresponds to a random index among possible indeces 
            const i = self.rand.uintLessThan(usize, indeces.items.len);
            tiles.set(indeces.items[i]);
        }

        fn getAdjacencies(self: *Self, p: GridIndex, k: usize) *const BitsetT {
            return &self.adjacencies[p][k];
        }

        /// =================
        /// = Debug methods
        /// =================

        fn getRandomCoords(self: *Self) GridIndex {
            return self.rand.uintLessThan(TileIndex, self.grid.len);
        }

        fn printAdjacencies(adjacencies: [][num_sides]BitsetT) void {
            for (adjacencies, 0..) |*arr, i| {
                std.debug.print("Adjacency for the {}th tile:\n", .{i});
                for (arr, 0..) |*b, k| {
                    std.debug.print("\t{}th side:\n", .{k});
                    printAdjacencyBitset(b, k);
                }
            }
        }

        fn printAdjacencyBitset(b: *BitsetT, k: usize) void {
            for (0..b.capacity()) |j| {
                if (b.isSet(j)) {
                    std.debug.print("\t\t{}th tile allowed on {}th side\n", .{j, (k+num_sides/2)%num_sides});
                }
            }
        }

        fn printBitset(b: *BitsetT) void {
            for (0..b.capacity()) |j| {
                if (b.isSet(j)) {
                    std.debug.print("\t\t{}th tile allowed\n", .{j});
                }
            }
        }

    };
}


const testing = std.testing;
const test_allocator = std.testing.allocator;

test "test toBuffer" {
    const t = SquareTile{.xpos = 0, .ypos = 1, .xneg = 2, .yneg = 3};
    const buf = Solver(SquareTile).toBuffer(t);

    for (0..4) |i| {
        try testing.expect(buf[i] == i);
    }
}