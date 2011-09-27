package kip.projects.cuda

import jcuda._
import jcuda.runtime._

class JCudaRuntimeTest extends App {
  val pointer = new Pointer()
  JCuda.cudaMalloc(pointer, 4)
  System.out.println("Pointer: "+pointer)
  JCuda.cudaFree(pointer)
}