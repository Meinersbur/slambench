
#include <stdint.h>

typedef unsigned int uint;
typedef unsigned short ushort;

typedef struct short2
{
	short x;
	short y;
} short2s;


/* Code to be extracted into pencil_kernel_core.cl. */
void initVolume_core(const uint x, const uint y, const uint z,
                            const uint v_size_x, const uint v_size_y, const uint v_size_z,
                            short2s *v_data,
                            const float dxVal, const float dyVal)
{
	short2s dVal;
	dVal.x = dxVal;
	dVal.y = dyVal;
	v_data[x + y * v_size_x + z * v_size_x * v_size_y] = dVal;
}

