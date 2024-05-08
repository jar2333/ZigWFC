const std = @import("std");

pub fn build(b: *std.Build) void {
    // Standard target options allows the person running `zig build` to choose
    // what target to build for. Here we do not override the defaults, which
    // means any target is allowed, and the default is native. Other options
    // for restricting supported target set are available.
    const target = b.standardTargetOptions(.{});

    // Options: Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall. Here we do not
    const optimize = b.standardOptimizeOption(.{});

    // ==============
    // = Zig library
    // ==============
    const lib = b.addStaticLibrary(.{
        .name = "wfc",
        .root_source_file = .{ .path = "src/wfc.zig" },
        .target = target,
        .optimize = optimize,
    });

    // This declares intent for the library to be installed into the standard
    // location when the user invokes the "install" step (the default step when
    // running `zig build`).
    b.installArtifact(lib);

    // ==============
    // = C library
    // ==============
    const c_lib = b.addStaticLibrary(.{
        .name = "wfc_c",
        .root_source_file = .{ .path = "src/wfc_c.zig" },
        .target = target,
        .optimize = optimize,
    });
    c_lib.linkLibC();
    b.installArtifact(c_lib);

    const demo = b.addExecutable(.{
        .name = "demo",
        .root_source_file = .{.path = "src/pipes.zig"},
        .target = target,
        .optimize = optimize,
    });

    b.installArtifact(demo);

    // This *creates* a Run step in the build graph, to be executed when another
    // step is evaluated that depends on it. The next line below will establish
    // such a dependency.
    const run_cmd = b.addRunArtifact(demo);

    // By making the run step depend on the install step, it will be run from the
    // installation directory rather than directly from within the cache directory.
    // This is not necessary, however, if the application depends on other installed
    // files, this ensures they will be present and in the expected location.
    run_cmd.step.dependOn(b.getInstallStep());

    // This allows the user to pass arguments to the application in the build
    // command itself, like this: `zig build run -- arg1 arg2 etc`
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    // This creates a build step. It will be visible in the `zig build --help` menu,
    // and can be selected like this: `zig build run`
    // This will evaluate the `run` step rather than the default, which is "install".
    const run_step = b.step("run", "Run the demo cli app");
    run_step.dependOn(&run_cmd.step);

//     // Creates a step for unit testing. This only builds the test executable
//     // but does not run it.
//     const main_tests = b.addTest(.{
//         .root_source_file = .{ .path = "src/main.zig" },
//         .target = target,
//         .optimize = optimize,
//     });

//     const run_main_tests = b.addRunArtifact(main_tests);

//     // This creates a build step. It will be visible in the `zig build --help` menu,
//     // and can be selected like this: `zig build test`
//     // This will evaluate the `test` step rather than the default, which is "install".
//     const test_step = b.step("test", "Run library tests");
//     test_step.dependOn(&run_main_tests.step);
}
