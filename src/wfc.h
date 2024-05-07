#ifndef ZIGWFC_H
#define ZIGWFC_H

#include <stdint.h>

struct wfc_SquareTile {
    uint8_t xpos;
    uint8_t ypos;
    uint8_t xneg;
    uint8_t yneg;
};

extern void *wfc_initSquareGridSolver(wfc_SquareTile *tiles, size_t tiles_size, uint64_t seed, int *err_code);
extern void wfc_freeSquareGridSolver(void *solver);

extern int wfc_solveSquareGrid(void *solver, uint8_t *grid, size_t grid_size, size_t width, size_t height);