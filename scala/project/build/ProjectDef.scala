import sbt._

class ProjectDef(info: ProjectInfo) extends DefaultProject(info) {
  /* val mvnScalaNLP = "ScalaNLP Maven2" at "http://repo.scalanlp.org/repo/" */
  /* val mvnscalaTools = "Scala-Tools" at "http://scala-tools.org/repo-snapshots/" */
  
  override def mainScalaSourcePath = "src" // default = src/main/scala
  override def mainJavaSourcePath = "src-java" // default = src/main/java
  
  override def libraryDependencies = Set(
    // "org.scalanlp" % "scalala" % "0.3.1-SNAPSHOT",
    // "jfree"  % "jfreechart"  % "1.0.12"
  ) ++ super.libraryDependencies
}
