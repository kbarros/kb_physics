import sbt._
import Keys._
import java.io.PrintWriter

// Key reference: http://harrah.github.com/xsbt/latest/sxr/Keys.scala.html

object HelloBuild extends Build {

  // ------------------------------
  // Launcher task
  //
  val Mklauncher = config("mklauncher") extend(Compile)
  val mklauncher = TaskKey[Unit]("mklauncher")
  val mklauncherTask = mklauncher <<= (target, fullClasspath in Runtime) map { (target, cp) =>
    def writeFile(file: File, str: String) {
      val writer = new PrintWriter(file)
      writer.println(str)
      writer.close()
    }
    val jarStrs = cp.map(_.data.toString).filter(_.endsWith(".jar"))
    val ensimeString = jarStrs.map("\""+_+"\"").mkString(":compile-jars("," ",")")
    writeFile((target / "ensime-classpath").asFile, ensimeString)
    val cpString = cp.map(_.data).mkString(":")
    val launchString = """
CLASSPATH="%s"
JOGL_LIBRARY_PATH="/Users/kbarros/dev/repo/scala/lib"
scala -J-Xmx1g -usejavacp -Djava.class.path="${CLASSPATH}" -Djava.library.path="${JOGL_LIBRARY_PATH}" "$@"
""".format(cpString)
    val targetFile = (target / "sbt-launcher").asFile
    writeFile(targetFile, launchString)
    targetFile.setExecutable(true)
  }
  
  // ------------------------------
  // Project definition
  //
  lazy val smatrix = RootProject(file("../../smatrix"))
  lazy val project = Project (
    "project",
    file ("."),
    settings = Defaults.defaultSettings ++ Seq(mklauncherTask)
  ) dependsOn smatrix
}
