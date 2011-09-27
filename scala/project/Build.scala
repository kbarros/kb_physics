import sbt._
import Keys._
import java.io.PrintWriter

// Key reference: http://harrah.github.com/xsbt/latest/sxr/Keys.scala.html

object HelloBuild extends Build {

  // ------------------------------
  // Launcher task
  //
  val mklauncher = TaskKey[Unit]("mklauncher")
  val mklauncherTask = mklauncher <<= (target, unmanagedBase, fullClasspath in Runtime) map { (target, ub, cp) =>
    val cpString  = cp map (_.data) mkString ":"
    val libString = ((ub ** "*") filter (_.isDirectory) get) mkString ":"
    val launchString = """
CLASSPATH="%s"
LIBRARY_PATH="%s"
scala -J-Xmx1g -usejavacp -Djava.class.path="${CLASSPATH}" -Djava.library.path="${LIBRARY_PATH}" "$@"
""".format(cpString, libString)
    val targetFile = (target / "sbt-launcher").asFile
    val writer = new PrintWriter(targetFile)
    writer.println(launchString)
    writer.close()
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
