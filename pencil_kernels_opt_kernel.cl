#include </home/meinersbur/src/SLAMBench_1_1/cl_kernel_vector.cl>
__kernel void kernel0(__global unsigned short *in, __global float *out, int outSize_x, int outSize_y, int inSize_x, int inSize_y, int ratio)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);
    int private_xr;
    int private_yr;

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    if (outSize_y >= 16 * b0 + t0 + 1) {
      private_xr = ((16 * b1 + t1) * (ratio));
      private_yr = ((16 * b0 + t0) * (ratio));
      out[(16 * b0 + t0) * outSize_x + (16 * b1 + t1)] = (in[private_yr * inSize_x + private_xr] / 1000.0f);
    }
}
__kernel void kernel1(float e_d, __global float *gaussian, __global float *in, __global float *out, int size_x, int size_y, int gaussianS, int r)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    if (in[(8 * b0 + t0) * size_x + (16 * b1 + t1)] == 0) {
      out[(8 * b0 + t0) * size_x + (16 * b1 + t1)] = 0;
    } else {
      out[(8 * b0 + t0) * size_x + (16 * b1 + t1)] = bilateralFilter_core((16 * b1 + t1), (8 * b0 + t0), (size_x), (size_y), (r), (gaussianS), e_d, in, gaussian);
    }
}
__kernel void kernel2(__global float *depth, __global const Matrix4 *invK, __global struct float3 *vertex, int imageSize_x, int imageSize_y)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    if (imageSize_y >= 16 * b0 + t0 + 1)
      vertex[(16 * b0 + t0) * imageSize_x + (16 * b1 + t1)] = depth2vertex_core((16 * b1 + t1), (16 * b0 + t0), (imageSize_x), (imageSize_y), depth, invK[0]);
}
__kernel void kernel3(__global struct float3 *in, __global struct float3 *out, int imageSize_x, int imageSize_y)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    out[(4 * b0 + t0) * imageSize_x + (4 * b1 + t1)] = vertex2normal_core((4 * b1 + t1), (4 * b0 + t0), (imageSize_x), (imageSize_y), in);
}
__kernel void kernel4(const float e_d, __global float *in, __global float *out, int outSize_x, int outSize_y, int inSize_x, int inSize_y, int r)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    out[(4 * b0 + t0) * outSize_x + (16 * b1 + t1)] = halfSampleRobustImage_core((16 * b1 + t1), (4 * b0 + t0), (outSize_x), (outSize_y), (inSize_x), (inSize_y), in, e_d, (r));
}
__kernel void kernel5(__global struct float3 *normal, __global struct uchar3 *out, int normalSize_x, int normalSize_y)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    out[(4 * b0 + t0) * normalSize_x + (16 * b1 + t1)] = renderNormal_core((16 * b1 + t1), (4 * b0 + t0), (normalSize_x), (normalSize_y), normal);
}
__kernel void kernel6(__global float *depth, const float farPlane, const float nearPlane, __global struct uchar4 *out, float rangeScale, int depthSize_x, int depthSize_y)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    out[(4 * b0 + t0) * depthSize_x + (16 * b1 + t1)] = renderDepth_core((16 * b1 + t1), (4 * b0 + t0), (depthSize_x), (depthSize_y), depth, nearPlane, farPlane, rangeScale);
}
__kernel void kernel7(__global struct TrackData *data, __global struct uchar4 *out, int outSize_x, int outSize_y)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    out[(8 * b0 + t0) * outSize_x + (32 * b1 + t1)] = renderTrack_core((32 * b1 + t1), (8 * b0 + t0), (outSize_x), (outSize_y), data);
}
__kernel void kernel8(__global const float3 *ambient, const float farPlane, const float largestep, __global const float3 *light, const float nearPlane, __global struct uchar4 *out, const float step, __global const Matrix4 *view, __global struct short2 *volume_data, __global const float3 *volume_dim, int volume_size_x, int volume_size_y, int volume_size_z, int depthSize_x, int depthSize_y)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    out[(8 * b0 + t0) * depthSize_x + (32 * b1 + t1)] = renderVolume_core((32 * b1 + t1), (8 * b0 + t0), (volume_size_x), (volume_size_y), (volume_size_z), volume_data, volume_dim[0], view[0], nearPlane, farPlane, step, largestep, light[0], ambient[0]);
}
__kernel void kernel9(__global const Matrix4 *Ttrack, const float dist_threshold, __global struct float3 *inNormal, __global struct float3 *inVertex, const float normal_threshold, __global struct TrackData *output, __global struct float3 *refNormal, __global struct float3 *refVertex, __global const Matrix4 *view, int refSize_x, int refSize_y, int inSize_x, int inSize_y)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    if (inSize_y >= 8 * b0 + t0 + 1 && inSize_x >= 32 * b1 + t1 + 1)
      output[(8 * b0 + t0) * refSize_x + (32 * b1 + t1)] = track_core((refSize_x), (refSize_y), output[(8 * b0 + t0) * refSize_x + (32 * b1 + t1)], inVertex[(8 * b0 + t0) * inSize_x + (32 * b1 + t1)], inNormal[(8 * b0 + t0) * inSize_x + (32 * b1 + t1)], refVertex, refNormal, Ttrack[0], view[0], dist_threshold, normal_threshold);
}
__kernel void kernel10(__global float *intrmdSums, int Jsize_x, int Jsize_y, int size_x, int size_y)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    for (int c4 = 0; c4 < size_x; c4 += 1)
      intrmdSums[(c4 * 8 + (4 * b0 + t0)) * 32 + (16 * b1 + t1)] = 0;
}
__kernel void kernel11(__global float *intrmdSums, __global struct TrackData *J, int Jsize_x, int Jsize_y, int size_x, int size_y)
{
    int b0 = get_group_id(0);
    int t0 = get_local_id(0);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    #define min(x,y)    ((x) < (y) ? (x) : (y))
    for (int c1 = 0; c1 < size_y - 4 * b0 - 3; c1 += 16) {
      for (int c3 = 4 * b0 + t0; c3 <= min(15, size_y - c1 - 1); c3 += 8)
        for (int c4 = 0; c4 < size_x; c4 += 1)
          reduce_core((intrmdSums + (c4 * 8 + (4 * b0 + t0)) * 32), J[(c1 + c3) * Jsize_x + c4]);
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    }
}
__kernel void kernel12(__global float *sums, int Jsize_x, int Jsize_y, int size_x, int size_y)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    sums[(4 * b0 + t0) * 32 + (16 * b1 + t1)] = 0;
}
__kernel void kernel13(__global float *intrmdSums, __global float *sums, int Jsize_x, int Jsize_y, int size_x, int size_y)
{
    int b0 = get_group_id(0), b1 = get_group_id(1);
    int t0 = get_local_id(0), t1 = get_local_id(1);

    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    for (int c4 = 0; c4 < size_x; c4 += 1)
      sums[t0 * 32 + (16 * b1 + t1)] += intrmdSums[(c4 * 8 + t0) * 32 + (16 * b1 + t1)];
}
