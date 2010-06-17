import sbt._

class ProjectDef(info: ProjectInfo) extends DefaultProject(info) with JoglProject {
  override def mainScalaSourcePath = "src" // default = src/main/scala
  override def mainJavaSourcePath = "src-java" // default = src/main/java
  
  // val jvmOptions = Seq("-Djava.library.path="+(managedDependencyPath / "compile"))
  // override def fork = forkRun(jvmOptions)
}
