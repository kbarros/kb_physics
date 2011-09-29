package kip.projects.cuda

import java.io._


object CuPreparePtx {

  def fromSrcString(cuCode: String): String = {
    import java.io._
    val cuFile = File.createTempFile("src", ".cu")
    val writer = new BufferedWriter(new FileWriter(cuFile))
    writer.write(cuCode);
    writer.close();
    fromSrcFile(cuFile.getPath)
  }

  /**
   * The extension of the given file name is replaced with "ptx".
   * If the file with the resulting name does not exist, it is
   * compiled from the given file using NVCC. The name of the
   * PTX file is returned.
   *
   * @param cuFileName The name of the .CU file
   * @return The name of the PTX file
   * @throws IOException If an I/O error occurs
   */
  def fromSrcFile(cuFileName: String): String = {
    var endIndex = cuFileName.lastIndexOf('.');
    if (endIndex == -1) {
      endIndex = cuFileName.length() - 1;
    }
    val ptxFileName = cuFileName.substring(0, endIndex + 1) + "ptx";
    val ptxFile = new File(ptxFileName);
    if (ptxFile.exists()) {
      return ptxFileName;
    }

    val cuFile = new File(cuFileName);
    if (!cuFile.exists()) {
      throw new IOException("Input file not found: " + cuFileName);
    }
    val modelString = "-m" + System.getProperty("sun.arch.data.model");
    val command =
      "nvcc " + modelString + " -ptx " +
        cuFile.getPath() + " -o " + ptxFileName;

    System.out.println("Executing\n" + command);
    val process = Runtime.getRuntime().exec(command);

    val errorMessage = new String(toByteArray(process.getErrorStream()));
    val outputMessage = new String(toByteArray(process.getInputStream()));
    var exitValue = 0;
    try {
      exitValue = process.waitFor();
    } catch {
      case e: InterruptedException =>
        Thread.currentThread().interrupt();
        throw new IOException(
          "Interrupted while waiting for nvcc output", e);
    }

    if (exitValue != 0) {
      System.out.println("nvcc process exitValue " + exitValue);
      System.out.println("errorMessage:\n" + errorMessage);
      System.out.println("outputMessage:\n" + outputMessage);
      throw new IOException(
        "Could not create .ptx file: " + errorMessage);
    }

    System.out.println("Finished creating PTX file");
    ptxFileName;
  }

  /**
   * Fully reads the given InputStream and returns it as a byte array
   *
   * @param inputStream The input stream to read
   * @return The byte array containing the data from the input stream
   * @throws IOException If an I/O error occurs
   */
  def toByteArray(inputStream: InputStream): Array[Byte] = {
    val baos = new ByteArrayOutputStream();
    val buffer = new Array[Byte](8192)
    var read = inputStream.read(buffer)
    while (read != -1) {
      baos.write(buffer, 0, read);
      read = inputStream.read(buffer)
    }
    baos.toByteArray();
  }
}