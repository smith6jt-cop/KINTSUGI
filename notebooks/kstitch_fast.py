import numpy as np
try:
    from numba import njit, prange
    _HAS_NUMBA = True
except Exception:
    _HAS_NUMBA = False


if _HAS_NUMBA:
    @njit(parallel=True, cache=True)
    def _paste_no_blend(stitched, tiles, y_pos, x_pos, tile_h, tile_w):
        for i in prange(tiles.shape[0]):
            y0 = int(y_pos[i]); x0 = int(x_pos[i])
            for yy in range(tile_h):
                for xx in range(tile_w):
                    stitched[y0+yy, x0+xx] = tiles[i, yy, xx]
        return stitched

    @njit(parallel=True, cache=True)
    def _paste_blend(accum, weight, tiles, y_pos, x_pos, tile_h, tile_w):
        for i in prange(tiles.shape[0]):
            y0 = int(y_pos[i]); x0 = int(x_pos[i])
            for yy in range(tile_h):
                for xx in range(tile_w):
                    accum[y0+yy, x0+xx] += tiles[i, yy, xx]
                    weight[y0+yy, x0+xx] += 1.0
        return accum, weight
else:
    def _paste_no_blend(stitched, tiles, y_pos, x_pos, tile_h, tile_w):
        for i in range(tiles.shape[0]):
            y0 = int(y_pos[i]); x0 = int(x_pos[i])
            stitched[y0:y0+tile_h, x0:x0+tile_w] = tiles[i]
        return stitched

    def _paste_blend(accum, weight, tiles, y_pos, x_pos, tile_h, tile_w):
        for i in range(tiles.shape[0]):
            y0 = int(y_pos[i]); x0 = int(x_pos[i])
            accum[y0:y0+tile_h, x0:x0+tile_w] += tiles[i]
            weight[y0:y0+tile_h, x0:x0+tile_w] += 1.0
        return accum, weight


def paste_no_blend(tiles, y_pos2, x_pos2, out_h, out_w, out_dtype):
    stitched = np.zeros((out_h, out_w), dtype=out_dtype)
    tile_h, tile_w = tiles.shape[1:3]
    _paste_no_blend(stitched, tiles, y_pos2, x_pos2, tile_h, tile_w)
    return stitched


def paste_blend_mean(tiles, y_pos2, x_pos2, out_h, out_w, out_dtype):
    tile_h, tile_w = tiles.shape[1:3]
    accum = np.zeros((out_h, out_w), dtype=np.float64)
    weight = np.zeros((out_h, out_w), dtype=np.float64)
    _paste_blend(accum, weight, tiles.astype(np.float64, copy=False), y_pos2, x_pos2, tile_h, tile_w)
    weight[weight == 0] = 1.0
    return (accum / weight).astype(out_dtype, copy=False)
