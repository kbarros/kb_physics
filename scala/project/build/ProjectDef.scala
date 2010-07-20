import sbt._

class ProjectDef(info: ProjectInfo) extends DefaultProject(info) with JoglProject {
  override def mainScalaSourcePath = "src" // default = src/main/scala
  override def mainJavaSourcePath = "src-java" // default = src/main/java
  
  val jvmOptions = Seq("-Djava.library.path="+(managedDependencyPath / "compile"))
  override def fork = forkRun(jvmOptions)
  
  lazy val mkrunner = task {
    val jlineJar = runClasspath.get.find(_.toString.contains("jline"))
    val toolClasspathStr = Path.makeString(buildScalaJars.get ++ jlineJar)
    val runClasspathStr  = Path.makeString(runClasspath.get)
    val scalaHomeStr = buildLibraryJar.asFile.getParentFile.getParent
    val scalaRunner =
"""[ -n "$JAVA_OPTS" ] || JAVA_OPTS="-Xmx256M -Xms32M"
for i
do
  case "$i" in
    -D*)
      JAVA_OPTS="$JAVA_OPTS $i" ;;
    *)
      ;;
  esac
done 
if [ -z "$JAVACMD" -a -n "$JAVA_HOME" -a -x "$JAVA_HOME/bin/java" ]; then
    JAVACMD="$JAVA_HOME/bin/java"
fi
TOOL_CLASSPATH="%s"
RUN_CLASSPATH="%s"
SCALA_HOME="%s"
exec "${JAVACMD:=java}" $JAVA_OPTS -cp "$TOOL_CLASSPATH" -Dscala.home="$SCALA_HOME" -Denv.emacs="$EMACS" scala.tools.nsc.MainGenericRunner -cp "$RUN_CLASSPATH" "$@"
""".format(toolClasspathStr, runClasspathStr, scalaHomeStr)
    val targetFile = ("target"/"bin"/"scala").asFile
    FileUtilities.write(targetFile, scalaRunner, log)
    targetFile.setExecutable(true)
    None
  }
  
  lazy val print = task {log.info("hello"); None}
}
