package kip.projects.cuda


import java.nio.Buffer
import java.nio.DoubleBuffer
import java.nio.FloatBuffer
import java.nio.IntBuffer
import jcuda.driver.JCudaDriver.cuCtxCreate
import jcuda.driver.JCudaDriver.cuDeviceComputeCapability
import jcuda.driver.JCudaDriver.cuDeviceGet
import jcuda.driver.JCudaDriver.cuDeviceGetCount
import jcuda.driver.JCudaDriver.cuDeviceGetName
import jcuda.driver.JCudaDriver.cuInit
import jcuda.driver.JCudaDriver.cuModuleGetFunction
import jcuda.driver.JCudaDriver.cuModuleLoad
import jcuda.driver.CUcontext
import jcuda.driver.CUdevice
import jcuda.driver.CUfunction
import jcuda.driver.CUmodule
import jcuda.jcusparse.JCusparse
import jcuda.runtime.JCuda.cudaMalloc
import jcuda.runtime.JCuda.cudaMemcpy
import jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost
import jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice
import jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToDevice
import jcuda.runtime.JCuda
import jcuda.Pointer
import jcuda.Sizeof


class JCudaWorld() {
  cuInit(0)
  JCusparse.setExceptionsEnabled(true);
  JCuda.setExceptionsEnabled(true);

  val device = new CUdevice();
  cuDeviceGet(device, 0);
  val context = new CUcontext();
  cuCtxCreate(context, 0, device);
    
  def destroy() {
  }
  
  var functions = Map[String, CUfunction]()
  
  def loadModule(ptxFileName: String, fnNames: Seq[String]) = {
    val module = new CUmodule();
    cuModuleLoad(module, ptxFileName);
    fnNames foreach { name =>
      val function = new CUfunction();
      cuModuleGetFunction(function, module, name);
      functions += name -> function
    }
  }
  
  def printDeviceProperties() {
    // Obtain the number of devices
    val deviceCountArray = Array(0)
    cuDeviceGetCount(deviceCountArray);
    val deviceCount = deviceCountArray(0)
    println("Found " + deviceCount + " devices");

    for (i <- 0 until deviceCount) {
      val device = new CUdevice();
      cuDeviceGet(device, i);

      // Obtain the device name
      val deviceName = new Array[Byte](1024)
      cuDeviceGetName(deviceName, deviceName.length, device);
      val name = deviceName.map(_.toChar).mkString

      // Obtain the compute capability
      val majorArray = Array(0);
      val minorArray = Array(0);
      cuDeviceComputeCapability(
          majorArray, minorArray, device);
      val major = majorArray(0);
      val minor = minorArray(0);

      println("Device " + i + ": " + name + " with Compute Capability " + major + "." + minor);
    }    
  }
  
  def storageBytes[T](x: T): Int = {
    x match {
      case xp: Array[Int]    => xp.size * Sizeof.INT
      case xp: Array[Float]  => xp.size * Sizeof.FLOAT
      case xp: Array[Double] => xp.size * Sizeof.DOUBLE
      case xp: IntBuffer     => xp.capacity * Sizeof.INT
      case xp: FloatBuffer   => xp.capacity * Sizeof.FLOAT
      case xp: DoubleBuffer  => xp.capacity * Sizeof.DOUBLE
    }
  }
  
  def pointerTo(x: Any): Pointer = {
    x match {
      case xp: Array[Int]    => Pointer.to(xp)
      case xp: Array[Float]  => Pointer.to(xp)
      case xp: Array[Double] => Pointer.to(xp)
      case xp: Buffer        => Pointer.to(xp)
      case xp: Pointer       => Pointer.to(xp)
      case _ => sys.error("Failed to create host pointer for object "+x)
    }
  }
  
  def allocDeviceArray(data_h: Any): Pointer = {
    val data_d = new Pointer()
    val bytes = storageBytes(data_h)
    cudaMalloc(data_d, bytes)
    cudaMemcpy(data_d, pointerTo(data_h), bytes, cudaMemcpyHostToDevice)
    data_d
  }
  
  def cpyHostToDevice(data_d: Pointer, data_h: Any) {
    cudaMemcpy(data_d, pointerTo(data_h), storageBytes(data_h), cudaMemcpyHostToDevice)
  }
  
  def cpyDeviceToHost(data_h: Any, data_d: Pointer) {
    cudaMemcpy(pointerTo(data_h), data_d, storageBytes(data_h), cudaMemcpyDeviceToHost)
  }
  
  def cpyDeviceToDevice(d1: Pointer, d2: Pointer, nbytes: Int) {
    cudaMemcpy(d1, d2, nbytes, cudaMemcpyDeviceToDevice)
  }
}