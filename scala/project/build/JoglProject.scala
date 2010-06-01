import sbt._

trait JoglProject extends BasicManagedProject {
  def jogl_vers = "1.1.1a"
  val Windows = "Windows(.*)".r
  
  def jogl_sig = (
    (system[String]("os.name").value, system[String]("os.arch").value) match {
      case ("Linux", "i386") => "linux" :: "i586" :: Nil
      case ("Mac OS X" | "Darwin", _) => "macosx" :: "universal" :: Nil
      case (Windows(_), "x86") => "windows" :: "i586" :: Nil
      case (name, arch) => name :: arch :: Nil
    }
  ) map { _.toLowerCase.split(" ").mkString } mkString "-"
  
  def jogl_loc = 
    "http://download.java.net/media/jogl/builds/archive/jsr-231-%s/jogl-%s-%s.zip" format
      (jogl_vers, jogl_vers, jogl_sig)
  
  val jogl = "net.java.dev" % ("jogl-" + jogl_sig) % jogl_vers % "provided->default" from jogl_loc
  
  override def updateAction = super.updateAction && task {
    import FileUtilities._
    val lib_compile = configurationPath(Configurations.Compile)
    val lib_provided = configurationPath(Configurations.Provided)
    val jogl_zip = outputPath / "jogl_zip"
    ((lib_provided * "jogl-*.zip").get.toList flatMap { file =>
      unzip(file, jogl_zip, "jogl-*/lib/*", log).left.toOption.orElse {
        FileUtilities.clean(file, log)
      }
    } match {
      case Nil => None
      case list => Some(list mkString "\n")
    }) orElse {
      val files = (jogl_zip ** "lib" ##) ** "*"
      println("Copying %s to %s " format (files, lib_compile))
      copy(files.get, lib_compile, log).left.toOption
    }
  }
}
