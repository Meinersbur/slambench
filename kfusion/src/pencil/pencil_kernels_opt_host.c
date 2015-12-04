#include <assert.h>
#include <stdio.h>
#include <prl_scop.h>
typedef unsigned int uint;
typedef unsigned short ushort;

/* OpenCL vector type definitions for PENCIL. */
struct float3 {
	float x;
	float y;
	float z;
};

struct float4 {
	float x;
	float y;
	float z;
	float w;
};

struct uchar3 {
	unsigned char x;
	unsigned char y;
	unsigned char z;
};

struct uchar4 {
	unsigned char x;
	unsigned char y;
	unsigned char z;
	unsigned char w;
};

struct short2 {
	short x;
	short y;
};

struct float2 {
	float x;
	float y;
};

struct uint2 {
	uint x;
	uint y;
};

struct Matrix4 {
	struct float4 data[4];
};

struct TrackData {
	int result;
	float error;
	float J[6];
};

typedef struct float3 float3;
typedef struct float4 float4;
typedef struct uchar3 uchar3;
typedef struct uchar4 uchar4;
typedef struct short2 short2;
typedef struct float2 float2;
typedef struct uint2 uint2;
typedef struct Matrix4 Matrix4;
typedef struct TrackData TrackData;

static float3 make_float3(float x, float y, float z) {
	float3 ret;
	ret.x = x;
	ret.y = y;
	ret.z = z;
	return ret;
}

static float3 c_rotate(const Matrix4 M, const float3 v)
{
	return make_float3(M.data[0].x * v.x + M.data[0].y * v.y + M.data[0].z * v.z,
	                   M.data[1].x * v.x + M.data[1].y * v.y + M.data[1].z * v.z,
	                   M.data[2].x * v.x + M.data[2].y * v.y + M.data[2].z * v.z);
}

void bilateralFilter_core_summary(int x, int y, int size_x, int size_y,
                                  int r, int gaussianS, float e_d,
                                  const float in[restrict const static size_y][size_x],
                                  const float gaussian[restrict const static gaussianS])
{
	const float center = in[y][x];
	for (int i = 0; i <= gaussianS; ++i) {
		const float factor = gaussian[i + r];
	}
}

float bilateralFilter_core(int x, int y, int size_x, int size_y,
                           int r, int gaussianS, float e_d,
                           const float in[restrict const static size_y][size_x],
                           const float gaussian[restrict const static gaussianS])
      __attribute__((pencil_access(bilateralFilter_core_summary)));

void integrateKernel_core_summary(const uint vol_size_x, const uint vol_size_y,
                                  const uint vol_size_z, const float3 vol_dim,
                                  short2 vol_data[restrict const static vol_size_z][vol_size_y][vol_size_x],
                                  const uint x, const uint y, uint depthSize_x, uint depthSize_y,
                                  const float depth[restrict const static depthSize_y][depthSize_x],
                                  const Matrix4 invTrack, const Matrix4 K,
                                  const float mu, const float maxweight,
                                  const float3 delta, const float3 cameraDelta)
{
	for (int z = 0; z <= vol_size_z; ++z) {
		const float depthVal = depth[y][x];
		const short2 volVal = vol_data[z][y][x];
		vol_data[z][y][x] = volVal;
	}
	for (int i = 0; i < depthSize_y; i++)
	{
		for (int j = 0; j < depthSize_x; ++j)
		{
			const float val = depth[i][j];
		}
	}
}

void integrateKernel_core(const uint vol_size_x, const uint vol_size_y,
                          const uint vol_size_z, const float3 vol_dim,
                          short2 vol_data[restrict const static vol_size_z][vol_size_y][vol_size_x],
                          const uint x, const uint y, uint depthSize_x, uint depthSize_y,
                          const float depth[restrict const static depthSize_y][depthSize_x],
                          const Matrix4 invTrack, const Matrix4 K,
                          const float mu, const float maxweight,
                          const float3 delta, const float3 cameraDelta)
     __attribute__((pencil_access(integrateKernel_core_summary)));

void initVolume_core_summary(const uint x, const uint y, const uint z,
                             const uint v_size_x, const uint v_size_y, const uint v_size_z,
                             short2 v_data[restrict const static v_size_z][v_size_y][v_size_x],
                             const float dxVal, const float dyVal)
{
	short2 temp;
	temp.x = dxVal;
	temp.y = dyVal;
	v_data[z][y][x] = temp;
}

void initVolume_core(const uint x, const uint y, const uint z,
                     const uint v_size_x, const uint v_size_y, const uint v_size_z,
                     short2 v_data[restrict const static v_size_z][v_size_y][v_size_x],
                     const float dxVal, const float dyVal)
     __attribute__((pencil_access(initVolume_core_summary)));
void depth2vertex_core_summary(uint x, uint y, uint imageSize_x, uint imageSize_y,
                               const float depth[restrict const static imageSize_y][imageSize_x],
                               const Matrix4 invK)
{
	const float depth_val = depth[y][x];
}

float3 depth2vertex_core(const uint x, const uint y,
                         const uint imageSize_x, const uint imageSize_y,
                         const float depth[restrict const static imageSize_y][imageSize_x],
                         const Matrix4 invK)
       __attribute__((pencil_access(depth2vertex_core_summary)));

void vertex2normal_core_summary(const uint x, const uint y,
                                const uint imageSize_x, const uint imageSize_y,
                                const float3 in[restrict const static imageSize_y][imageSize_x])
{
	const float3 left = in[y][x];
}
float3 vertex2normal_core(const uint x, const uint y,
                          const uint imageSize_x, const uint imageSize_y,
                          const float3 in[restrict const static imageSize_y][imageSize_x])
       __attribute__((pencil_access(vertex2normal_core_summary)));

void halfSampleRobustImage_core_summary(const uint x, const uint y,
                                        const uint outSize_x, const uint outSize_y,
                                        const uint inSize_x, const uint inSize_y,
                                        const float in[restrict const static inSize_y][inSize_x],
                                        const float e_d, const int r)
{
	const float center = in[2*y][2*x];
}

float halfSampleRobustImage_core(const uint x, const uint y,
                                 const uint outSize_x, const uint outSize_y,
                                 const uint inSize_x, const uint inSize_y,
                                 const float in[restrict const static inSize_y][inSize_x],
                                 const float e_d, const int r)
      __attribute__((pencil_access(halfSampleRobustImage_core_summary)));

void renderNormal_core_summary(const uint x, const uint y,
                               const uint normalSize_x, const uint normalSize_y,
                               const float3 normal[restrict const static normalSize_y][normalSize_x])
{
	const float3 n = normal[y][x];
}

uchar3 renderNormal_core(const uint x, const uint y,
                         const uint normalSize_x, const uint normalSize_y,
                         const float3 normal[restrict const static normalSize_y][normalSize_x])
       __attribute__((pencil_access(renderNormal_core_summary)));

void renderDepth_core_summary(const uint x, const uint y,
                              const uint depthSize_x, const uint depthSize_y,
                              const float depth[restrict const static depthSize_y][depthSize_x],
                              const float nearPlane, const float farPlane,
                              const float rangeScale)
{
	const float d = depth[y][x];
}

uchar4 renderDepth_core(const uint x, const uint y,
                        const uint depthSize_x, const uint depthSize_y,
                        const float depth[restrict const static depthSize_y][depthSize_x],
                        const float nearPlane, const float farPlane,
                        const float rangeScale)
       __attribute__((pencil_access(renderDepth_core_summary)));

void renderTrack_core_summary(const uint x, const uint y,
                              const uint outSize_x, const uint outSize_y,
                              const TrackData data[restrict const static outSize_y][outSize_x])
{
	int test = data[y][x].result;
}

uchar4 renderTrack_core(const uint x, const uint y,
                        const uint outSize_x, const uint outSize_y,
                        const TrackData data[restrict const static outSize_y][outSize_x])
       __attribute__((pencil_access(renderTrack_core_summary)));

void renderVolume_core_summary(const uint x, const uint y,
                               const uint volume_size_x, const uint volume_size_y, const uint volume_size_z,
                               const short2 volume_data[restrict const static volume_size_z][volume_size_y][volume_size_x],
                               const float3 volume_dim, const Matrix4 view,
                               const float nearPlane, const float farPlane,
                               const float step, const float largestep,
                               const float3 light, const float3 ambient)
{
	const short2 d = volume_data[x+y][y][x];
}

uchar4 renderVolume_core(const uint x, const uint y,
                         const uint volume_size_x, const uint volume_size_y, const uint volume_size_z,
                         const short2 volume_data[restrict const static volume_size_z][volume_size_y][volume_size_x],
                         const float3 volume_dim, const Matrix4 view,
                         const float nearPlane, const float farPlane,
                         const float step, const float largestep,
                         const float3 light, const float3 ambient)
       __attribute__((pencil_access(renderVolume_core_summary)));

void raycast_core_summary(const uint x, const uint y,
                          const uint inputSize_x, const uint inputSize_y,
                          float3 vertex[restrict const static inputSize_y][inputSize_x],
                          float3 normal[restrict const static inputSize_y][inputSize_x],
                          const uint integration_size_x, const uint integration_size_y, const uint integration_size_z,
                          const short2 integration_data[restrict const static integration_size_z][integration_size_y][integration_size_x],
                          const float3 integration_dim, const Matrix4 view,
                          const float nearPlane, const float farPlane,
                          const float step, const float largestep)
{
	const float3 test;
	const float3 norm;
	const short2 integrVal = integration_data[x+y][y][x];
	vertex[y][x] = test;
	normal[y][x] = norm;
}

void raycast_core(const uint x, const uint y,
                  const uint inputSize_x, const uint inputSize_y,
                  float3 vertex[restrict const static inputSize_y][inputSize_x],
                  float3 normal[restrict const static inputSize_y][inputSize_x],
                  const uint integration_size_x, const uint integration_size_y, const uint integration_size_z,
                  const short2 integration_data[restrict const static integration_size_z][integration_size_y][integration_size_x],
                  const float3 integration_dim, const Matrix4 view,
                  const float nearPlane, const float farPlane,
                  const float step, const float largestep)
     __attribute__((pencil_access(raycast_core_summary)));

void track_core_summary(uint refSize_x, uint refSize_y, const TrackData output,
                        const float3 inVertex, const float3 inNormal,
                        const float3 refVertex[restrict const static refSize_y][refSize_x],
                        const float3 refNormal[restrict const static refSize_y][refSize_x],
                        const Matrix4 Ttrack, const Matrix4 view,
                        const float dist_threshold,
                        const float normal_threshold)
{
	const uint refx;
	const uint refy;
	const float3 vertVal = refVertex[refy][refx];
	const float3 normVal = refNormal[refy][refx];
}

TrackData track_core(uint refSize_x, uint refSize_y, const TrackData output,
                     const float3 inVertex, const float3 inNormal,
                     const float3 refVertex[restrict const static refSize_y][refSize_x],
                     const float3 refNormal[restrict const static refSize_y][refSize_x],
                     const Matrix4 Ttrack, const Matrix4 view,
                     const float dist_threshold,
                     const float normal_threshold)
          __attribute__((pencil_access(track_core_summary)));

void reduce_core_summary(float sums[restrict const static 32], TrackData row)
{
	for (int z = 0; z < 32; ++z) {
		const float adjustVal = row.J[z/6];
		const float tempVal = sums[z];
		sums[z] = adjustVal + tempVal;
	}
}
void reduce_core(float sums[restrict const static 32], TrackData row)
     __attribute__((pencil_access(reduce_core_summary)));


int mm2meters_pencil(uint outSize_x, uint outSize_y,
					 float out[restrict const static outSize_y][outSize_x],
					 uint inSize_x, uint inSize_y,
					 const ushort in[restrict const static inSize_y][inSize_x],
					 int ratio)
{
	#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	{
	  #define max(x,y)    ((x) > (y) ? (x) : (y))
	  static prl_mem __ppcg_mem_in;
	  static prl_mem __ppcg_mem_out;
	  
	  static prl_program __ppcg_program;
	  static prl_scop __ppcg_scop;
	  prl_scop_instance __ppcg_scopinst;
	  
	  __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	  prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	  __ppcg_mem_in = prl_scop_get_mem(__ppcg_scopinst, in, max(sizeof(unsigned short), (inSize_y) * (inSize_x) * sizeof(unsigned short)), "in");
	  __ppcg_mem_out = prl_scop_get_mem(__ppcg_scopinst, out, max(sizeof(float), (outSize_y) * (outSize_x) * sizeof(float)), "out");
	  
	  {
	    prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_in);
	    {
	      size_t grid_size[2] = {floord(outSize_y + 15, 16), outSize_x / 16};
	      size_t block_size[2] = {16, 16};
	      static prl_kernel __ppcg_kernel0;
	      prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel0, __ppcg_program, "kernel0");
	      struct prl_kernel_call_arg __ppcg_kernel0_args[] = {
	        { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_in },
	        { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_out },
	        { .type = prl_kernel_call_arg_value, .data = &outSize_x, .size = sizeof(outSize_x) }, 
	        { .type = prl_kernel_call_arg_value, .data = &outSize_y, .size = sizeof(outSize_y) }, 
	        { .type = prl_kernel_call_arg_value, .data = &inSize_x, .size = sizeof(inSize_x) }, 
	        { .type = prl_kernel_call_arg_value, .data = &inSize_y, .size = sizeof(inSize_y) }, 
	        { .type = prl_kernel_call_arg_value, .data = &ratio, .size = sizeof(ratio) }, 
	       };
	       prl_scop_call(__ppcg_scopinst, __ppcg_kernel0, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel0_args)/sizeof(__ppcg_kernel0_args[0]), __ppcg_kernel0_args);
	     }
	     
	     prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_out);
	   }
	   prl_scop_leave(__ppcg_scopinst);
	 }
	return 0;
}

int bilateralFilter_pencil(int size_x, int size_y,
						   float out[restrict const static size_y][size_x],
						   const float in[restrict const static size_y][size_x],
						   uint2 size, int gaussianS,
						   const float gaussian[restrict const static gaussianS],
						   float e_d, int r)
{
	 #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	 {
	   #define max(x,y)    ((x) > (y) ? (x) : (y))
	   static prl_mem __ppcg_mem_gaussian;
	   static prl_mem __ppcg_mem_in;
	   static prl_mem __ppcg_mem_out;
	   
	   static prl_program __ppcg_program;
	   static prl_scop __ppcg_scop;
	   prl_scop_instance __ppcg_scopinst;
	   
	   __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	   prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	   __ppcg_mem_gaussian = prl_scop_get_mem(__ppcg_scopinst, gaussian, max(sizeof(float), (gaussianS) * sizeof(float)), "gaussian");
	   __ppcg_mem_in = prl_scop_get_mem(__ppcg_scopinst, in, max(sizeof(float), (size_y) * (size_x) * sizeof(float)), "in");
	   __ppcg_mem_out = prl_scop_get_mem(__ppcg_scopinst, out, max(sizeof(float), (size_y) * (size_x) * sizeof(float)), "out");
	   
	   {
	     if (gaussianS >= r + 1)
	       prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_gaussian);
	     prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_in);
	     prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_out);
	     {
	       size_t grid_size[2] = {size_y / 8, size_x / 16};
	       size_t block_size[2] = {8, 16};
	       static prl_kernel __ppcg_kernel1;
	       prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel1, __ppcg_program, "kernel1");
	       struct prl_kernel_call_arg __ppcg_kernel1_args[] = {
	         { .type = prl_kernel_call_arg_value, .data = &e_d, .size = sizeof(e_d) }, 
	         { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_gaussian },
	         { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_in },
	         { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_out },
	         { .type = prl_kernel_call_arg_value, .data = &size_x, .size = sizeof(size_x) }, 
	         { .type = prl_kernel_call_arg_value, .data = &size_y, .size = sizeof(size_y) }, 
	         { .type = prl_kernel_call_arg_value, .data = &gaussianS, .size = sizeof(gaussianS) }, 
	         { .type = prl_kernel_call_arg_value, .data = &r, .size = sizeof(r) }, 
	        };
	        prl_scop_call(__ppcg_scopinst, __ppcg_kernel1, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel1_args)/sizeof(__ppcg_kernel1_args[0]), __ppcg_kernel1_args);
	      }
	      
	      prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_out);
	    }
	    prl_scop_leave(__ppcg_scopinst);
	  }
	return 0;
}

int initVolume_pencil(const uint v_size_x, const uint v_size_y, const uint v_size_z,
					  short2 v_data[restrict const static v_size_z][v_size_y][v_size_x],
					  const float2 d)
{
	  /* ppcg generated CPU code */
	  
	  #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	  for (int c0 = 0; c0 < v_size_x; c0 += 1)
	    for (int c1 = 0; c1 < v_size_y; c1 += 1)
	      for (int c2 = 0; c2 < v_size_z; c2 += 1)
	        initVolume_core((c0), (c1), (c2), (v_size_x), (v_size_y), (v_size_z), v_data, d.x * 32766.0f, d.y);
	return 0;
}

int integrateKernel_pencil(const uint vol_size_x, const uint vol_size_y,
						   const uint vol_size_z, const float3 vol_dim,
						   short2 vol_data[restrict const static vol_size_z][vol_size_y][vol_size_x],
						   uint depthSize_x, uint depthSize_y,
						   const float depth[restrict const static depthSize_y][depthSize_x],
						   const Matrix4 invTrack, const Matrix4 K,
						   const float mu, const float maxweight)
{
	const float3 delta = c_rotate(invTrack,
								  make_float3(0, 0, vol_dim.z / vol_size_z));
	const float3 cameraDelta = c_rotate(K, delta);
	  /* ppcg generated CPU code */
	  
	  #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	  for (int c0 = 0; c0 < vol_size_y; c0 += 1)
	    for (int c1 = 0; c1 < vol_size_x; c1 += 1)
	      integrateKernel_core((vol_size_x), (vol_size_y), (vol_size_z), vol_dim, vol_data, (c1), (c0), (depthSize_x), (depthSize_y), depth, invTrack, K, mu, maxweight, delta, cameraDelta);
	return 0;
}

int depth2vertex_pencil(uint imageSize_x, uint imageSize_y,
						float3 vertex[restrict const static imageSize_y][imageSize_x],
						const float depth[restrict const static imageSize_y][imageSize_x],
						const Matrix4 invK)
{
	  #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	  {
	    #define max(x,y)    ((x) > (y) ? (x) : (y))
	    static prl_mem __ppcg_mem_depth;
	    static prl_mem __ppcg_mem_invK;
	    static prl_mem __ppcg_mem_vertex;
	    
	    static prl_program __ppcg_program;
	    static prl_scop __ppcg_scop;
	    prl_scop_instance __ppcg_scopinst;
	    
	    __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	    prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	    __ppcg_mem_depth = prl_scop_get_mem(__ppcg_scopinst, depth, max(sizeof(float), (imageSize_y) * (imageSize_x) * sizeof(float)), "depth");
	    __ppcg_mem_invK = prl_scop_get_mem(__ppcg_scopinst, &invK, sizeof(const Matrix4), "invK");
	    __ppcg_mem_vertex = prl_scop_get_mem(__ppcg_scopinst, vertex, max(sizeof(struct float3), (imageSize_y) * (imageSize_x) * sizeof(struct float3)), "vertex");
	    
	    {
	      prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_depth);
	      prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_invK);
	      {
	        size_t grid_size[2] = {floord(imageSize_y + 15, 16), imageSize_x / 16};
	        size_t block_size[2] = {16, 16};
	        static prl_kernel __ppcg_kernel2;
	        prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel2, __ppcg_program, "kernel2");
	        struct prl_kernel_call_arg __ppcg_kernel2_args[] = {
	          { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_depth },
	          { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_invK },
	          { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_vertex },
	          { .type = prl_kernel_call_arg_value, .data = &imageSize_x, .size = sizeof(imageSize_x) }, 
	          { .type = prl_kernel_call_arg_value, .data = &imageSize_y, .size = sizeof(imageSize_y) }, 
	         };
	         prl_scop_call(__ppcg_scopinst, __ppcg_kernel2, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel2_args)/sizeof(__ppcg_kernel2_args[0]), __ppcg_kernel2_args);
	       }
	       
	       prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_vertex);
	     }
	     prl_scop_leave(__ppcg_scopinst);
	   }
	return 0;
}

int vertex2normal_pencil(uint imageSize_x, uint imageSize_y,
						 float3 out[restrict const static imageSize_y][imageSize_x],
						 const float3 in[restrict const static imageSize_y][imageSize_x])
{
	   #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	   {
	     #define max(x,y)    ((x) > (y) ? (x) : (y))
	     static prl_mem __ppcg_mem_in;
	     static prl_mem __ppcg_mem_out;
	     
	     static prl_program __ppcg_program;
	     static prl_scop __ppcg_scop;
	     prl_scop_instance __ppcg_scopinst;
	     
	     __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	     prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	     __ppcg_mem_in = prl_scop_get_mem(__ppcg_scopinst, in, max(sizeof(struct float3), (imageSize_y) * (imageSize_x) * sizeof(struct float3)), "in");
	     __ppcg_mem_out = prl_scop_get_mem(__ppcg_scopinst, out, max(sizeof(struct float3), (imageSize_y) * (imageSize_x) * sizeof(struct float3)), "out");
	     
	     {
	       prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_in);
	       {
	         size_t grid_size[2] = {imageSize_y / 4, imageSize_x / 4};
	         size_t block_size[2] = {4, 4};
	         static prl_kernel __ppcg_kernel3;
	         prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel3, __ppcg_program, "kernel3");
	         struct prl_kernel_call_arg __ppcg_kernel3_args[] = {
	           { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_in },
	           { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_out },
	           { .type = prl_kernel_call_arg_value, .data = &imageSize_x, .size = sizeof(imageSize_x) }, 
	           { .type = prl_kernel_call_arg_value, .data = &imageSize_y, .size = sizeof(imageSize_y) }, 
	          };
	          prl_scop_call(__ppcg_scopinst, __ppcg_kernel3, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel3_args)/sizeof(__ppcg_kernel3_args[0]), __ppcg_kernel3_args);
	        }
	        
	        prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_out);
	      }
	      prl_scop_leave(__ppcg_scopinst);
	    }
	return 0;
}

int halfSampleRobustImage_pencil(uint outSize_x, uint outSize_y,
								 uint inSize_x, uint inSize_y,
								 float out[restrict const static outSize_y][outSize_x],
								 const float in[restrict const static inSize_y][inSize_x],
								 const float e_d, const int r)
{
	    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	    {
	      #define max(x,y)    ((x) > (y) ? (x) : (y))
	      static prl_mem __ppcg_mem_in;
	      static prl_mem __ppcg_mem_out;
	      
	      static prl_program __ppcg_program;
	      static prl_scop __ppcg_scop;
	      prl_scop_instance __ppcg_scopinst;
	      
	      __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	      prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	      __ppcg_mem_in = prl_scop_get_mem(__ppcg_scopinst, in, max(sizeof(float), (inSize_y >= 1 && 4294967294 * ((4294967295 * inSize_y + 4294967297) / 4294967296) >= 4294967293 * inSize_y + 2 * floord(4294967294 * outSize_y + inSize_y + 1, 4294967296) && (inSize_y + 4294967294) % 4294967296 <= 4294967294 ? 2 * floord(4294967294 * outSize_y + inSize_y + 1, 4294967296) - 1 : 2 * floord(inSize_y + 1, 2) - 1) * (inSize_x) * sizeof(float)), "in");
	      __ppcg_mem_out = prl_scop_get_mem(__ppcg_scopinst, out, max(sizeof(float), (outSize_y) * (outSize_x) * sizeof(float)), "out");
	      
	      {
	        if ((inSize_x >= 1 && 4294967296 * ((4294967293 * inSize_y + 2 * ((4294967294 * outSize_y + inSize_y + 1) / 4294967296) + 4294967293) / 4294967294) >= 4294967295 * inSize_y + 4294967297) || (inSize_x >= 1 && inSize_y >= 1 && 4294967294 * ((4294967295 * inSize_y + 4294967297) / 4294967296) >= 4294967293 * inSize_y + 2 * ((4294967294 * outSize_y + inSize_y + 1) / 4294967296) && (inSize_y + 4294967294) % 4294967296 <= 4294967294))
	          prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_in);
	        {
	          size_t grid_size[2] = {outSize_y / 4, outSize_x / 16};
	          size_t block_size[2] = {4, 16};
	          static prl_kernel __ppcg_kernel4;
	          prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel4, __ppcg_program, "kernel4");
	          struct prl_kernel_call_arg __ppcg_kernel4_args[] = {
	            { .type = prl_kernel_call_arg_value, .data = &e_d, .size = sizeof(e_d) }, 
	            { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_in },
	            { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_out },
	            { .type = prl_kernel_call_arg_value, .data = &outSize_x, .size = sizeof(outSize_x) }, 
	            { .type = prl_kernel_call_arg_value, .data = &outSize_y, .size = sizeof(outSize_y) }, 
	            { .type = prl_kernel_call_arg_value, .data = &inSize_x, .size = sizeof(inSize_x) }, 
	            { .type = prl_kernel_call_arg_value, .data = &inSize_y, .size = sizeof(inSize_y) }, 
	            { .type = prl_kernel_call_arg_value, .data = &r, .size = sizeof(r) }, 
	           };
	           prl_scop_call(__ppcg_scopinst, __ppcg_kernel4, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel4_args)/sizeof(__ppcg_kernel4_args[0]), __ppcg_kernel4_args);
	         }
	         
	         prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_out);
	       }
	       prl_scop_leave(__ppcg_scopinst);
	     }
	return 0;
}

int renderNormal_pencil(uint normalSize_x, uint normalSize_y,
						uchar3 out[restrict const static normalSize_y][normalSize_x],
						const float3 normal[restrict const static normalSize_y][normalSize_x])
{
	     #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	     {
	       #define max(x,y)    ((x) > (y) ? (x) : (y))
	       static prl_mem __ppcg_mem_normal;
	       static prl_mem __ppcg_mem_out;
	       
	       static prl_program __ppcg_program;
	       static prl_scop __ppcg_scop;
	       prl_scop_instance __ppcg_scopinst;
	       
	       __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	       prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	       __ppcg_mem_normal = prl_scop_get_mem(__ppcg_scopinst, normal, max(sizeof(struct float3), (normalSize_y) * (normalSize_x) * sizeof(struct float3)), "normal");
	       __ppcg_mem_out = prl_scop_get_mem(__ppcg_scopinst, out, max(sizeof(struct uchar3), (normalSize_y) * (normalSize_x) * sizeof(struct uchar3)), "out");
	       
	       {
	         prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_normal);
	         {
	           size_t grid_size[2] = {normalSize_y / 4, normalSize_x / 16};
	           size_t block_size[2] = {4, 16};
	           static prl_kernel __ppcg_kernel5;
	           prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel5, __ppcg_program, "kernel5");
	           struct prl_kernel_call_arg __ppcg_kernel5_args[] = {
	             { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_normal },
	             { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_out },
	             { .type = prl_kernel_call_arg_value, .data = &normalSize_x, .size = sizeof(normalSize_x) }, 
	             { .type = prl_kernel_call_arg_value, .data = &normalSize_y, .size = sizeof(normalSize_y) }, 
	            };
	            prl_scop_call(__ppcg_scopinst, __ppcg_kernel5, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel5_args)/sizeof(__ppcg_kernel5_args[0]), __ppcg_kernel5_args);
	          }
	          
	          prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_out);
	        }
	        prl_scop_leave(__ppcg_scopinst);
	      }
	return 0;
}

int renderDepth_pencil(uint depthSize_x, uint depthSize_y,
					   uchar4 out[restrict const static depthSize_y][depthSize_x],
					   const float depth[restrict const static depthSize_y][depthSize_x],
					   const float nearPlane, const float farPlane)
{
	float rangeScale = 1 / (farPlane - nearPlane);
	      #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	      {
	        #define max(x,y)    ((x) > (y) ? (x) : (y))
	        static prl_mem __ppcg_mem_depth;
	        static prl_mem __ppcg_mem_out;
	        
	        static prl_program __ppcg_program;
	        static prl_scop __ppcg_scop;
	        prl_scop_instance __ppcg_scopinst;
	        
	        __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	        prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	        __ppcg_mem_depth = prl_scop_get_mem(__ppcg_scopinst, depth, max(sizeof(float), (depthSize_y) * (depthSize_x) * sizeof(float)), "depth");
	        __ppcg_mem_out = prl_scop_get_mem(__ppcg_scopinst, out, max(sizeof(struct uchar4), (depthSize_y) * (depthSize_x) * sizeof(struct uchar4)), "out");
	        
	        {
	          prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_depth);
	          {
	            size_t grid_size[2] = {depthSize_y / 4, depthSize_x / 16};
	            size_t block_size[2] = {4, 16};
	            static prl_kernel __ppcg_kernel6;
	            prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel6, __ppcg_program, "kernel6");
	            struct prl_kernel_call_arg __ppcg_kernel6_args[] = {
	              { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_depth },
	              { .type = prl_kernel_call_arg_value, .data = &farPlane, .size = sizeof(farPlane) }, 
	              { .type = prl_kernel_call_arg_value, .data = &nearPlane, .size = sizeof(nearPlane) }, 
	              { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_out },
	              { .type = prl_kernel_call_arg_value, .data = &rangeScale, .size = sizeof(rangeScale) }, 
	              { .type = prl_kernel_call_arg_value, .data = &depthSize_x, .size = sizeof(depthSize_x) }, 
	              { .type = prl_kernel_call_arg_value, .data = &depthSize_y, .size = sizeof(depthSize_y) }, 
	             };
	             prl_scop_call(__ppcg_scopinst, __ppcg_kernel6, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel6_args)/sizeof(__ppcg_kernel6_args[0]), __ppcg_kernel6_args);
	           }
	           
	           prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_out);
	         }
	         prl_scop_leave(__ppcg_scopinst);
	       }
	return 0;
}

int renderTrack_pencil(uint outSize_x, uint outSize_y,
					   uchar4 out[restrict const static outSize_y][outSize_x],
					   const TrackData data[restrict const static outSize_y][outSize_x])
{
	       #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	       {
	         #define max(x,y)    ((x) > (y) ? (x) : (y))
	         static prl_mem __ppcg_mem_data;
	         static prl_mem __ppcg_mem_out;
	         
	         static prl_program __ppcg_program;
	         static prl_scop __ppcg_scop;
	         prl_scop_instance __ppcg_scopinst;
	         
	         __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	         prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	         __ppcg_mem_data = prl_scop_get_mem(__ppcg_scopinst, data, max(sizeof(struct TrackData), (outSize_y) * (outSize_x) * sizeof(struct TrackData)), "data");
	         __ppcg_mem_out = prl_scop_get_mem(__ppcg_scopinst, out, max(sizeof(struct uchar4), (outSize_y) * (outSize_x) * sizeof(struct uchar4)), "out");
	         
	         {
	           prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_data);
	           {
	             size_t grid_size[2] = {outSize_y / 8, outSize_x / 32};
	             size_t block_size[2] = {8, 32};
	             static prl_kernel __ppcg_kernel7;
	             prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel7, __ppcg_program, "kernel7");
	             struct prl_kernel_call_arg __ppcg_kernel7_args[] = {
	               { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_data },
	               { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_out },
	               { .type = prl_kernel_call_arg_value, .data = &outSize_x, .size = sizeof(outSize_x) }, 
	               { .type = prl_kernel_call_arg_value, .data = &outSize_y, .size = sizeof(outSize_y) }, 
	              };
	              prl_scop_call(__ppcg_scopinst, __ppcg_kernel7, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel7_args)/sizeof(__ppcg_kernel7_args[0]), __ppcg_kernel7_args);
	            }
	            
	            prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_out);
	          }
	          prl_scop_leave(__ppcg_scopinst);
	        }
	return 0;
}

int renderVolume_pencil(uint depthSize_x, uint depthSize_y,
						uchar4 out[restrict const static depthSize_y][depthSize_x],
						const uint volume_size_x, const uint volume_size_y, const uint volume_size_z,
						const short2 volume_data[restrict const static volume_size_z][volume_size_y][volume_size_x],
						const float3 volume_dim, const Matrix4 view,
						const float nearPlane, const float farPlane,
						const float step, const float largestep,
						const float3 light, const float3 ambient)
{
	        #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	        {
	          #define max(x,y)    ((x) > (y) ? (x) : (y))
	          static prl_mem __ppcg_mem_ambient;
	          static prl_mem __ppcg_mem_light;
	          static prl_mem __ppcg_mem_out;
	          static prl_mem __ppcg_mem_view;
	          static prl_mem __ppcg_mem_volume_data;
	          static prl_mem __ppcg_mem_volume_dim;
	          
	          static prl_program __ppcg_program;
	          static prl_scop __ppcg_scop;
	          prl_scop_instance __ppcg_scopinst;
	          
	          __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	          prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	          __ppcg_mem_ambient = prl_scop_get_mem(__ppcg_scopinst, &ambient, sizeof(const float3), "ambient");
	          __ppcg_mem_light = prl_scop_get_mem(__ppcg_scopinst, &light, sizeof(const float3), "light");
	          __ppcg_mem_out = prl_scop_get_mem(__ppcg_scopinst, out, max(sizeof(struct uchar4), (depthSize_y) * (depthSize_x) * sizeof(struct uchar4)), "out");
	          __ppcg_mem_view = prl_scop_get_mem(__ppcg_scopinst, &view, sizeof(const Matrix4), "view");
	          __ppcg_mem_volume_data = prl_scop_get_mem(__ppcg_scopinst, volume_data, max(sizeof(struct short2), ((volume_size_y >= 1 && volume_size_z + 4294967296 >= volume_size_y + depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && 4294967295 * depthSize_y + 18446744073709551616 * floord(volume_size_z + 4294967295, 4294967296) >= 4294967296 * volume_size_z + 18446744069414584320 && 4294967296 * floord(volume_size_x + 4294967295 * depthSize_x, 4294967296) >= 4294967295 * depthSize_x + 32) || (volume_size_y >= 1 && depthSize_x + depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967297 && volume_size_z + 4294967296 >= volume_size_y + depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && 4294967296 * volume_size_z + 18446744069414584312 >= 4294967295 * depthSize_y + 18446744073709551616 * floord(volume_size_z + 4294967295, 4294967296) && volume_size_x >= ((volume_size_x + 4294967295 * depthSize_x) % 4294967296) + 32) || (volume_size_y >= 1 && volume_size_z >= 1 && depthSize_y >= volume_size_y + 1 && ((volume_size_z - 1) % 4294967296) + 1 >= depthSize_x + depthSize_y && 4294967296 * floord(volume_size_x + 4294967295 * depthSize_x, 4294967296) >= 4294967295 * depthSize_x + 32) ? volume_size_y + depthSize_x - 1 : (depthSize_y + 4294967296 * floord(volume_size_y + 4294967295, 4294967296) >= volume_size_y + 4294967296 && depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967296 && volume_size_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967296 * floord(volume_size_y + 4294967295, 4294967296) + 1) || (volume_size_y >= 1 && (volume_size_y - 1) % 4294967296 >= depthSize_y && depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967297 && depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967296) || (volume_size_z + 4294967295 >= depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) + 4294967296 * floord(volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295, 4294967296) >= volume_size_y + 4294967296 * volume_size_z + 8589934592 && volume_size_y + 4294967295 * volume_size_z + 4294967295 >= 4294967296 * floord(volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295, 4294967296)) || (depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967297 && volume_size_z + 4294967295 >= depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && 4294967296 * floord(volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295, 4294967296) >= 4294967295 * volume_size_z + 4294967297 && volume_size_y + 4294967296 * volume_size_z + 8589934591 >= depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) + 4294967296 * floord(volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295, 4294967296)) || (depthSize_y + 4294967296 * floord(volume_size_y + 4294967295, 4294967296) >= volume_size_y + 4294967296 && depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967296 && volume_size_z + 4294967296 * floord(volume_size_y + 4294967295, 4294967296) >= volume_size_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && 4294967296 * floord(volume_size_x + volume_size_y + 4294967295 * volume_size_z + 4294967295, 4294967296) >= volume_size_y + 4294967295 * volume_size_z + 4294967296) || (volume_size_z + 4294967295 >= depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && 4294967296 * floord(volume_size_x + volume_size_y + 4294967295 * volume_size_z + 4294967295, 4294967296) >= volume_size_y + 4294967295 * volume_size_z + 4294967296 && depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) + 4294967296 * floord(volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295, 4294967296) >= volume_size_y + 4294967296 * volume_size_z + 8589934592 && 4294967296 * floord(volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295, 4294967296) >= volume_size_y + 4294967295 * volume_size_z + 4294967296) || (depthSize_x + depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967297 && volume_size_z + 4294967295 >= depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && volume_size_z + 4294967296 >= depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) + 4294967296 * floord(volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295, 4294967296) >= 4294967296 * volume_size_z + 8589934600 && volume_size_y + 4294967296 * volume_size_z + 8589934591 >= depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) + 4294967296 * floord(volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295, 4294967296) && 4294967296 * floord(volume_size_x + 4294967295 * volume_size_z + depthSize_y + 4294967295, 4294967296) >= 4294967295 * volume_size_z + depthSize_y + 4294967296) || (volume_size_y + 4294967295 >= depthSize_y + 4294967296 * floord(volume_size_y + 4294967295, 4294967296) && depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967296 && volume_size_z + 4294967296 >= depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && 4294967296 * floord(volume_size_x + 4294967295 * volume_size_z + depthSize_y + 4294967295, 4294967296) >= 4294967295 * volume_size_z + depthSize_y + 4294967296) ? volume_size_z : (volume_size_y >= 1 && depthSize_x >= volume_size_x + 1 && volume_size_z + 4294967296 >= volume_size_y + depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && 4294967295 * depthSize_y + 18446744073709551616 * floord(volume_size_z + 4294967295, 4294967296) >= 4294967296 * volume_size_z + 18446744069414584320) || (volume_size_y >= 1 && depthSize_x >= volume_size_x + 1 && depthSize_x + depthSize_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967297 && volume_size_z + 4294967296 >= volume_size_y + depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) && 4294967296 * volume_size_z + 18446744069414584312 >= 4294967295 * depthSize_y + 18446744073709551616 * floord(volume_size_z + 4294967295, 4294967296)) || (volume_size_y >= 1 && depthSize_y >= volume_size_y && volume_size_y + depthSize_x + 4294967296 * floord(volume_size_z + 4294967295, 4294967296) >= volume_size_z + 4294967297 && volume_size_z + 4294967296 >= volume_size_x + volume_size_y + 4294967296 * floord(volume_size_z + 4294967295, 4294967296)) || (volume_size_y >= 1 && volume_size_z >= 1 && depthSize_x >= volume_size_x + 1 && depthSize_y >= volume_size_y + 1 && ((volume_size_z - 1) % 4294967296) + 1 >= depthSize_x + depthSize_y) ? volume_size_x + volume_size_y - 1 : volume_size_z >= 1 && ((volume_size_z - 1) % 4294967296) + 1 >= depthSize_x + depthSize_y && 4294967296 * floord(volume_size_y + 4294967295 * depthSize_y, 4294967296) >= 4294967295 * depthSize_y + 8 && 4294967296 * floord(volume_size_x + 4294967295 * depthSize_x, 4294967296) >= 4294967295 * depthSize_x + 32 ? depthSize_x + depthSize_y - 1 : volume_size_x + depthSize_y - 1) * (volume_size_y) * (volume_size_x) * sizeof(struct short2)), "volume_data");
	          __ppcg_mem_volume_dim = prl_scop_get_mem(__ppcg_scopinst, &volume_dim, sizeof(const float3), "volume_dim");
	          
	          {
	            prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_ambient);
	            prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_light);
	            prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_view);
	            if ((volume_size_x >= 1 && volume_size_y >= 1 && depthSize_x >= volume_size_x + 1 && ((volume_size_z + 4294967295) % 4294967296) + 1 >= volume_size_y + depthSize_x && 4294967295 * depthSize_y >= 4294967296 * ((volume_size_z + 4294967295) % 4294967296)) || (volume_size_x >= 1 && volume_size_y >= 1 && depthSize_y >= volume_size_y && volume_size_y + depthSize_x >= ((volume_size_z + 4294967295) % 4294967296) + 2 && ((volume_size_z + 4294967295) % 4294967296) + 1 >= volume_size_x + volume_size_y) || (volume_size_y >= 1 && (volume_size_y - 1) % 4294967296 >= depthSize_y && depthSize_x >= ((volume_size_z + 4294967295) % 4294967296) + 1 && ((volume_size_z + 4294967295) % 4294967296) + 1 >= depthSize_y && volume_size_x >= ((volume_size_x + 4294967295 * volume_size_z + depthSize_y + 4294967295) % 4294967296) + 1) || (volume_size_x >= 1 && depthSize_y >= ((volume_size_y + 4294967295) % 4294967296) + 1 && depthSize_x >= ((volume_size_z + 4294967295) % 4294967296) + 1 && (volume_size_y + 4294967295) % 4294967296 >= ((volume_size_z + 4294967295) % 4294967296) + 1) || (volume_size_x >= 1 && volume_size_y >= 1 && depthSize_x >= volume_size_x + 1 && depthSize_x + depthSize_y >= ((volume_size_z + 4294967295) % 4294967296) + 2 && ((volume_size_z + 4294967295) % 4294967296) + 1 >= volume_size_y + depthSize_x && 4294967296 * ((volume_size_z + 4294967295) % 4294967296) >= 4294967295 * depthSize_y + 8) || (volume_size_x >= 1 && volume_size_y >= 1 && volume_size_z >= 1 && depthSize_x >= volume_size_x + 1 && depthSize_y >= volume_size_y + 1 && ((volume_size_z - 1) % 4294967296) + 1 >= depthSize_x + depthSize_y) || (volume_size_x >= 1 && volume_size_y >= 1 && (volume_size_y - 1) % 4294967296 >= depthSize_y && depthSize_y >= ((volume_size_z + 4294967295) % 4294967296) + 2 && depthSize_x >= ((volume_size_z + 4294967295) % 4294967296) + 1) || (volume_size_y >= 1 && ((volume_size_z + 4294967295) % 4294967296) + 1 >= volume_size_y + depthSize_x && 4294967295 * depthSize_y >= 4294967296 * ((volume_size_z + 4294967295) % 4294967296) && volume_size_x >= ((volume_size_x + 4294967295 * depthSize_x) % 4294967296) + 32) || (volume_size_y >= 1 && depthSize_x + depthSize_y >= ((volume_size_z + 4294967295) % 4294967296) + 2 && ((volume_size_z + 4294967295) % 4294967296) + 1 >= volume_size_y + depthSize_x && 4294967296 * ((volume_size_z + 4294967295) % 4294967296) >= 4294967295 * depthSize_y + 8 && volume_size_x >= ((volume_size_x + 4294967295 * depthSize_x) % 4294967296) + 32) || (volume_size_x >= 1 && (volume_size_z + 4294967295) % 4294967296 >= depthSize_x && depthSize_x + depthSize_y >= ((volume_size_z + 4294967295) % 4294967296) + ((volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295) % 4294967296) + 2 && (volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295) % 4294967296 >= depthSize_x) || (volume_size_x >= 1 && depthSize_x + depthSize_y >= ((volume_size_z + 4294967295) % 4294967296) + 2 && (volume_size_z + 4294967295) % 4294967296 >= depthSize_x && ((volume_size_z + 4294967295) % 4294967296) + 1 >= volume_size_x + depthSize_y && volume_size_y + depthSize_x + depthSize_y >= ((volume_size_z + 4294967295) % 4294967296) + ((volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295) % 4294967296) + 10 && ((volume_size_z + 4294967295) % 4294967296) + ((volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295) % 4294967296) + 1 >= depthSize_x + depthSize_y) || (depthSize_y >= ((volume_size_y + 4294967295) % 4294967296) + 1 && depthSize_x >= ((volume_size_z + 4294967295) % 4294967296) + 1 && (volume_size_z + 4294967295) % 4294967296 >= (volume_size_y + 4294967295) % 4294967296 && volume_size_x >= ((volume_size_x + volume_size_y + 4294967295 * volume_size_z + 4294967295) % 4294967296) + 1) || ((volume_size_z + 4294967295) % 4294967296 >= depthSize_x && volume_size_x >= ((volume_size_x + volume_size_y + 4294967295 * volume_size_z + 4294967295) % 4294967296) + 1 && depthSize_x + depthSize_y >= ((volume_size_z + 4294967295) % 4294967296) + ((volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295) % 4294967296) + 2 && depthSize_x >= ((volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295) % 4294967296) + 1) || (volume_size_y >= 1 && volume_size_z >= 1 && depthSize_y >= volume_size_y + 1 && ((volume_size_z - 1) % 4294967296) + 1 >= depthSize_x + depthSize_y && volume_size_x >= ((volume_size_x + 4294967295 * depthSize_x) % 4294967296) + 32) || (volume_size_x >= 1 && volume_size_z >= 1 && depthSize_x >= volume_size_x + 1 && ((volume_size_z - 1) % 4294967296) + 1 >= depthSize_x + depthSize_y && volume_size_y >= ((volume_size_y + 4294967295 * depthSize_y) % 4294967296) + 8) || (volume_size_x >= 1 && volume_size_y >= 1 && (volume_size_y - 1) % 4294967296 >= depthSize_y && depthSize_x >= ((volume_size_z + 4294967295) % 4294967296) + 1 && ((volume_size_z + 4294967295) % 4294967296) + 1 >= volume_size_x + depthSize_y) || (volume_size_x >= 1 && depthSize_y >= ((volume_size_z + 4294967295) % 4294967296) + 2 && (volume_size_z + 4294967295) % 4294967296 >= depthSize_x && volume_size_y + depthSize_x >= ((volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295) % 4294967296) + 2 && ((volume_size_z + 4294967295) % 4294967296) + ((volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295) % 4294967296) + 1 >= depthSize_x + depthSize_y) || (volume_size_z >= 1 && ((volume_size_z - 1) % 4294967296) + 1 >= depthSize_x + depthSize_y && volume_size_x >= ((volume_size_x + 4294967295 * depthSize_x) % 4294967296) + 32 && volume_size_y >= ((volume_size_y + 4294967295 * depthSize_y) % 4294967296) + 8) || (depthSize_x + depthSize_y >= ((volume_size_z + 4294967295) % 4294967296) + 2 && (volume_size_z + 4294967295) % 4294967296 >= depthSize_x && ((volume_size_z + 4294967295) % 4294967296) + 1 >= depthSize_y && volume_size_y + depthSize_x + depthSize_y >= ((volume_size_z + 4294967295) % 4294967296) + ((volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295) % 4294967296) + 10 && ((volume_size_z + 4294967295) % 4294967296) + ((volume_size_y + 4294967295 * volume_size_z + depthSize_x + 4294967295) % 4294967296) + 1 >= depthSize_x + depthSize_y && volume_size_x >= ((volume_size_x + 4294967295 * volume_size_z + depthSize_y + 4294967295) % 4294967296) + 1))
	              prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_volume_data);
	            prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_volume_dim);
	            {
	              size_t grid_size[2] = {depthSize_y / 8, depthSize_x / 32};
	              size_t block_size[2] = {8, 32};
	              static prl_kernel __ppcg_kernel8;
	              prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel8, __ppcg_program, "kernel8");
	              struct prl_kernel_call_arg __ppcg_kernel8_args[] = {
	                { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_ambient },
	                { .type = prl_kernel_call_arg_value, .data = &farPlane, .size = sizeof(farPlane) }, 
	                { .type = prl_kernel_call_arg_value, .data = &largestep, .size = sizeof(largestep) }, 
	                { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_light },
	                { .type = prl_kernel_call_arg_value, .data = &nearPlane, .size = sizeof(nearPlane) }, 
	                { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_out },
	                { .type = prl_kernel_call_arg_value, .data = &step, .size = sizeof(step) }, 
	                { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_view },
	                { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_volume_data },
	                { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_volume_dim },
	                { .type = prl_kernel_call_arg_value, .data = &volume_size_x, .size = sizeof(volume_size_x) }, 
	                { .type = prl_kernel_call_arg_value, .data = &volume_size_y, .size = sizeof(volume_size_y) }, 
	                { .type = prl_kernel_call_arg_value, .data = &volume_size_z, .size = sizeof(volume_size_z) }, 
	                { .type = prl_kernel_call_arg_value, .data = &depthSize_x, .size = sizeof(depthSize_x) }, 
	                { .type = prl_kernel_call_arg_value, .data = &depthSize_y, .size = sizeof(depthSize_y) }, 
	               };
	               prl_scop_call(__ppcg_scopinst, __ppcg_kernel8, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel8_args)/sizeof(__ppcg_kernel8_args[0]), __ppcg_kernel8_args);
	             }
	             
	             prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_out);
	           }
	           prl_scop_leave(__ppcg_scopinst);
	         }
	return 0;
}

int raycast_pencil(uint inputSize_x, uint inputSize_y,
				   float3 vertex[restrict const static inputSize_y][inputSize_x],
				   float3 normal[restrict const static inputSize_y][inputSize_x],
				   const uint integration_size_x, const uint integration_size_y, const uint integration_size_z,
				   const short2 integration_data[restrict const static integration_size_z][integration_size_y][integration_size_x],
				   const float3 integration_dim, const Matrix4 view,
				   const float nearPlane, const float farPlane,
				   const float step, const float largestep)
{
	         /* ppcg generated CPU code */
	         
	         #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	         for (int c0 = 0; c0 < inputSize_y; c0 += 1)
	           for (int c1 = 0; c1 < inputSize_x; c1 += 1)
	             raycast_core((c1), (c0), (inputSize_x), (inputSize_y), vertex, normal, (integration_size_x), (integration_size_y), (integration_size_z), integration_data, integration_dim, view, nearPlane, farPlane, step, largestep);
	return 0;
}

int track_pencil(uint refSize_x, uint refSize_y, uint inSize_x, uint inSize_y,
				 TrackData output[restrict const static refSize_y][refSize_x],
				 const float3 inVertex[restrict const static inSize_y][inSize_x],
				 const float3 inNormal[restrict const static inSize_y][inSize_x],
				 const float3 refVertex[restrict const static refSize_y][refSize_x],
				 const float3 refNormal[restrict const static refSize_y][refSize_x],
				 const Matrix4 Ttrack, const Matrix4 view,
				 const float dist_threshold, const float normal_threshold)
{
	         #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	         {
	           #define max(x,y)    ((x) > (y) ? (x) : (y))
	           static prl_mem __ppcg_mem_Ttrack;
	           static prl_mem __ppcg_mem_inNormal;
	           static prl_mem __ppcg_mem_inVertex;
	           static prl_mem __ppcg_mem_output;
	           static prl_mem __ppcg_mem_refNormal;
	           static prl_mem __ppcg_mem_refVertex;
	           static prl_mem __ppcg_mem_view;
	           
	           static prl_program __ppcg_program;
	           static prl_scop __ppcg_scop;
	           prl_scop_instance __ppcg_scopinst;
	           
	           __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	           prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	           __ppcg_mem_Ttrack = prl_scop_get_mem(__ppcg_scopinst, &Ttrack, sizeof(const Matrix4), "Ttrack");
	           __ppcg_mem_inNormal = prl_scop_get_mem(__ppcg_scopinst, inNormal, max(sizeof(struct float3), (inSize_y) * (inSize_x) * sizeof(struct float3)), "inNormal");
	           __ppcg_mem_inVertex = prl_scop_get_mem(__ppcg_scopinst, inVertex, max(sizeof(struct float3), (inSize_y) * (inSize_x) * sizeof(struct float3)), "inVertex");
	           __ppcg_mem_output = prl_scop_get_mem(__ppcg_scopinst, output, max(sizeof(struct TrackData), (refSize_y >= inSize_y ? inSize_y : refSize_y) * (refSize_x) * sizeof(struct TrackData)), "output");
	           __ppcg_mem_refNormal = prl_scop_get_mem(__ppcg_scopinst, refNormal, max(sizeof(struct float3), (refSize_y) * (refSize_x) * sizeof(struct float3)), "refNormal");
	           __ppcg_mem_refVertex = prl_scop_get_mem(__ppcg_scopinst, refVertex, max(sizeof(struct float3), (refSize_y) * (refSize_x) * sizeof(struct float3)), "refVertex");
	           __ppcg_mem_view = prl_scop_get_mem(__ppcg_scopinst, &view, sizeof(const Matrix4), "view");
	           
	           {
	             prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_Ttrack);
	             prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_inNormal);
	             prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_inVertex);
	             if ((refSize_x >= 1 && refSize_y >= inSize_y) || (refSize_x >= 1 && refSize_y >= 1 && inSize_y >= refSize_y + 1)) {
	               prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_output);
	               if (refSize_y >= 1) {
	                 prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_refNormal);
	                 prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_refVertex);
	               }
	             }
	             prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_view);
	             {
	               size_t grid_size[2] = {floord(inSize_y + 7, 8), floord(inSize_x + 31, 32)};
	               size_t block_size[2] = {8, 32};
	               static prl_kernel __ppcg_kernel9;
	               prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel9, __ppcg_program, "kernel9");
	               struct prl_kernel_call_arg __ppcg_kernel9_args[] = {
	                 { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_Ttrack },
	                 { .type = prl_kernel_call_arg_value, .data = &dist_threshold, .size = sizeof(dist_threshold) }, 
	                 { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_inNormal },
	                 { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_inVertex },
	                 { .type = prl_kernel_call_arg_value, .data = &normal_threshold, .size = sizeof(normal_threshold) }, 
	                 { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_output },
	                 { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_refNormal },
	                 { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_refVertex },
	                 { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_view },
	                 { .type = prl_kernel_call_arg_value, .data = &refSize_x, .size = sizeof(refSize_x) }, 
	                 { .type = prl_kernel_call_arg_value, .data = &refSize_y, .size = sizeof(refSize_y) }, 
	                 { .type = prl_kernel_call_arg_value, .data = &inSize_x, .size = sizeof(inSize_x) }, 
	                 { .type = prl_kernel_call_arg_value, .data = &inSize_y, .size = sizeof(inSize_y) }, 
	                };
	                prl_scop_call(__ppcg_scopinst, __ppcg_kernel9, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel9_args)/sizeof(__ppcg_kernel9_args[0]), __ppcg_kernel9_args);
	              }
	              
	              if ((refSize_x >= 1 && refSize_y >= inSize_y) || (refSize_x >= 1 && refSize_y >= 1 && inSize_y >= refSize_y + 1))
	                prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_output);
	            }
	            prl_scop_leave(__ppcg_scopinst);
	          }
	return 0;
}

int reduce_pencil(float sums[restrict const static 8][32], const uint Jsize_x, const uint Jsize_y,
				  TrackData J[restrict const static Jsize_y][Jsize_x],
				  const uint size_x, const uint size_y)
{
	          #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
	          {
	            #define max(x,y)    ((x) > (y) ? (x) : (y))
	            static prl_mem __ppcg_mem_intrmdSums;
	            static prl_mem __ppcg_mem_J;
	            static prl_mem __ppcg_mem_sums;
	            
	            static prl_program __ppcg_program;
	            static prl_scop __ppcg_scop;
	            prl_scop_instance __ppcg_scopinst;
	            
	            __ppcg_scopinst = prl_scop_enter(&__ppcg_scop);
	            prl_scop_program_from_file(__ppcg_scopinst, &__ppcg_program, "pencil_kernels_opt_kernel.cl", "-I.");

	            __ppcg_mem_intrmdSums = prl_scop_get_mem(__ppcg_scopinst, NULL, max(sizeof(float), (size_x) * (8) * (32) * sizeof(float)), "intrmdSums");
	            __ppcg_mem_J = prl_scop_get_mem(__ppcg_scopinst, J, max(sizeof(struct TrackData), (Jsize_y >= size_y + 60 ? size_y : Jsize_y) * (Jsize_x) * sizeof(struct TrackData)), "J");
	            __ppcg_mem_sums = prl_scop_get_mem(__ppcg_scopinst, sums, (8) * (32) * sizeof(float), "sums");
	            
	            {
	              if ((Jsize_x >= 160 && Jsize_y >= size_y + 60) || (Jsize_x >= 160 && Jsize_y >= 120 && size_y >= Jsize_y))
	                prl_scop_host_to_device(__ppcg_scopinst, __ppcg_mem_J);
	              {
	                size_t grid_size[2] = {2, 2};
	                size_t block_size[2] = {4, 16};
	                static prl_kernel __ppcg_kernel10;
	                prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel10, __ppcg_program, "kernel10");
	                struct prl_kernel_call_arg __ppcg_kernel10_args[] = {
	                  { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_intrmdSums },
	                  { .type = prl_kernel_call_arg_value, .data = &Jsize_x, .size = sizeof(Jsize_x) }, 
	                  { .type = prl_kernel_call_arg_value, .data = &Jsize_y, .size = sizeof(Jsize_y) }, 
	                  { .type = prl_kernel_call_arg_value, .data = &size_x, .size = sizeof(size_x) }, 
	                  { .type = prl_kernel_call_arg_value, .data = &size_y, .size = sizeof(size_y) }, 
	                 };
	                 prl_scop_call(__ppcg_scopinst, __ppcg_kernel10, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel10_args)/sizeof(__ppcg_kernel10_args[0]), __ppcg_kernel10_args);
	               }
	               
	               {
	                 size_t grid_size[1] = {2};
	                 size_t block_size[1] = {4};
	                 static prl_kernel __ppcg_kernel11;
	                 prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel11, __ppcg_program, "kernel11");
	                 struct prl_kernel_call_arg __ppcg_kernel11_args[] = {
	                   { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_intrmdSums },
	                   { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_J },
	                   { .type = prl_kernel_call_arg_value, .data = &Jsize_x, .size = sizeof(Jsize_x) }, 
	                   { .type = prl_kernel_call_arg_value, .data = &Jsize_y, .size = sizeof(Jsize_y) }, 
	                   { .type = prl_kernel_call_arg_value, .data = &size_x, .size = sizeof(size_x) }, 
	                   { .type = prl_kernel_call_arg_value, .data = &size_y, .size = sizeof(size_y) }, 
	                  };
	                  prl_scop_call(__ppcg_scopinst, __ppcg_kernel11, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel11_args)/sizeof(__ppcg_kernel11_args[0]), __ppcg_kernel11_args);
	                }
	                
	                {
	                  size_t grid_size[2] = {2, 2};
	                  size_t block_size[2] = {4, 16};
	                  static prl_kernel __ppcg_kernel12;
	                  prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel12, __ppcg_program, "kernel12");
	                  struct prl_kernel_call_arg __ppcg_kernel12_args[] = {
	                    { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_sums },
	                    { .type = prl_kernel_call_arg_value, .data = &Jsize_x, .size = sizeof(Jsize_x) }, 
	                    { .type = prl_kernel_call_arg_value, .data = &Jsize_y, .size = sizeof(Jsize_y) }, 
	                    { .type = prl_kernel_call_arg_value, .data = &size_x, .size = sizeof(size_x) }, 
	                    { .type = prl_kernel_call_arg_value, .data = &size_y, .size = sizeof(size_y) }, 
	                   };
	                   prl_scop_call(__ppcg_scopinst, __ppcg_kernel12, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel12_args)/sizeof(__ppcg_kernel12_args[0]), __ppcg_kernel12_args);
	                 }
	                 
	                 {
	                   size_t grid_size[2] = {1, 2};
	                   size_t block_size[2] = {8, 16};
	                   static prl_kernel __ppcg_kernel13;
	                   prl_scop_init_kernel(__ppcg_scopinst, &__ppcg_kernel13, __ppcg_program, "kernel13");
	                   struct prl_kernel_call_arg __ppcg_kernel13_args[] = {
	                     { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_intrmdSums },
	                     { .type = prl_kernel_call_arg_mem, .mem = __ppcg_mem_sums },
	                     { .type = prl_kernel_call_arg_value, .data = &Jsize_x, .size = sizeof(Jsize_x) }, 
	                     { .type = prl_kernel_call_arg_value, .data = &Jsize_y, .size = sizeof(Jsize_y) }, 
	                     { .type = prl_kernel_call_arg_value, .data = &size_x, .size = sizeof(size_x) }, 
	                     { .type = prl_kernel_call_arg_value, .data = &size_y, .size = sizeof(size_y) }, 
	                    };
	                    prl_scop_call(__ppcg_scopinst, __ppcg_kernel13, sizeof(grid_size)/sizeof(grid_size[0]), grid_size,  sizeof(block_size)/sizeof(block_size[0]), block_size, sizeof(__ppcg_kernel13_args)/sizeof(__ppcg_kernel13_args[0]), __ppcg_kernel13_args);
	                  }
	                  
	                  prl_scop_device_to_host(__ppcg_scopinst, __ppcg_mem_sums);
	                }
	                prl_scop_leave(__ppcg_scopinst);
	              }
	return 0;
}
