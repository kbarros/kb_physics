package kip.projects.cuda

import jcuda.driver.JCudaDriver._;

import java.io._;

import jcuda._;
import jcuda.driver._

/**
 * This is a sample class demonstrating how to use the JCuda driver
 * bindings to load and execute a CUDA vector addition kernel.
 * The sample reads a CUDA file, compiles it to a PTX file
 * using NVCC, loads the PTX file as a module and executes
 * the kernel function. <br />
 */
object JCudaVectorAdd {
  /**
   * Entry point of this sample
   *
   * @param args Not used
   * @throws IOException If an IO error occurs
   */
  def main(args: Array[String]) {
    // Enable exceptions and omit all subsequent error checks
    JCudaDriver.setExceptionsEnabled(true);

    // Create the PTX file by calling the NVCC
//    val ptxFileName = preparePtxFile("JCudaVectorAddKernel.cu");
    val ptxFileName = CuPreparePtx.fromSrcString("""
extern "C"
__global__ void add(int n, float *a, float *b, float *sum)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i<n)
    {
        sum[i] = a[i] + b[i];
    }

}        
"""
    )

    // Initialize the driver and create a context for the first device.
    cuInit(0);
    val device = new CUdevice();
    cuDeviceGet(device, 0);
    val context = new CUcontext();
    cuCtxCreate(context, 0, device);

    // Load the ptx file.
    val module = new CUmodule();
    cuModuleLoad(module, ptxFileName);

    // Obtain a function pointer to the "add" function.
    val function = new CUfunction();
    cuModuleGetFunction(function, module, "add");

    val numElements = 100000;

    // Allocate and fill the host input data
    val hostInputA = new Array[Float](numElements)
    val hostInputB = new Array[Float](numElements)
    for (i <- 0 until numElements) {
      hostInputA(i) = i;
      hostInputB(i) = i;
    }

    // Allocate the device input data, and copy the
    // host input data to the device
    val deviceInputA = new CUdeviceptr();
    cuMemAlloc(deviceInputA, numElements * Sizeof.FLOAT);
    cuMemcpyHtoD(deviceInputA, Pointer.to(hostInputA),
      numElements * Sizeof.FLOAT);
    val deviceInputB = new CUdeviceptr();
    cuMemAlloc(deviceInputB, numElements * Sizeof.FLOAT);
    cuMemcpyHtoD(deviceInputB, Pointer.to(hostInputB),
      numElements * Sizeof.FLOAT);

    // Allocate device output memory
    val deviceOutput = new CUdeviceptr();
    cuMemAlloc(deviceOutput, numElements * Sizeof.FLOAT);

    // Set up the kernel parameters: A pointer to an array
    // of pointers which point to the actual values.
    val kernelParameters = Pointer.to(
      Pointer.to(Array(numElements)),
      Pointer.to(deviceInputA),
      Pointer.to(deviceInputB),
      Pointer.to(deviceOutput));

    // Call the kernel function.
    val blockSizeX = 256;
    val gridSizeX = Math.ceil(numElements.toDouble / blockSizeX).toInt;
    cuLaunchKernel(function,
      gridSizeX, 1, 1, // Grid dimension
      blockSizeX, 1, 1, // Block dimension
      0, null, // Shared memory size and stream
      kernelParameters, null // Kernel- and extra parameters
      );
    cuCtxSynchronize();

    // Allocate host output memory and copy the device output
    // to the host.
    val hostOutput = new Array[Float](numElements)
    cuMemcpyDtoH(Pointer.to(hostOutput), deviceOutput,
      numElements * Sizeof.FLOAT);

    // Verify the result
    var passed = true;
    for (i <- 0 until numElements; if passed) {
      val expected = i + i;
      if (Math.abs(hostOutput(i) - expected) > 1e-5) {
        println(
          "At index " + i + " found " + hostOutput(i) +
            " but expected " + expected);
        passed = false;
      }
    }
    println("Test " + (if (passed) "PASSED" else "FAILED"));

    // Clean up.
    cuMemFree(deviceInputA);
    cuMemFree(deviceInputB);
    cuMemFree(deviceOutput);
  }
}
