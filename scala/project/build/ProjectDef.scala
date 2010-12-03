import sbt._

class ProjectDef(info: ProjectInfo) extends DefaultProject(info) with JoglProject with AkkaProject {
  val AkkaRepo = "Akka Repository" at "http://scalablesolutions.se/akka/repository"
  
  val json = "com.twitter" % "json" % "2.1.3"
  
  override def compileOptions: Seq[CompileOption] = Deprecation :: Unchecked :: Nil
  
  // non-standard source directories
  override def mainScalaSourcePath = "src" // default = src/main/scala
  override def mainJavaSourcePath = "src-java" // default = src/main/java

  // needed for Jogl dependencies
  val javaLibraryPath = Path.makeString(Seq(managedDependencyPath / "compile"))
  val jvmOptions = Seq("-Djava.library.path="+javaLibraryPath)
  override def fork = forkRun(jvmOptions)
  
  // creates a runner script "sbt-scala" with project automatically on classpath
  lazy val mkrunner = task {
    val jlineJar = runClasspath.get.find(_.toString.contains("jline"))
    val toolClasspathStr = Path.makeString(buildScalaJars.get ++ jlineJar)
    val runClasspathStr  = Path.makeString(runClasspath.get ++ buildScalaJars.get)
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
JAVA_LIBRARY_PATH="%s"
exec "${JAVACMD:=java}" $JAVA_OPTS -cp "$TOOL_CLASSPATH" -Dscala.home="$SCALA_HOME" -Dscala.class.path="$RUN_CLASSPATH" -Denv.emacs="$EMACS" -Djava.library.path="$JAVA_LIBRARY_PATH" scala.tools.nsc.MainGenericRunner -cp "$RUN_CLASSPATH" "$@"
""".format(toolClasspathStr, runClasspathStr, scalaHomeStr, javaLibraryPath)
    val targetFile = ("target"/"bin"/"sbt-scala").asFile
    FileUtilities.write(targetFile, scalaRunner, log)
    targetFile.setExecutable(true)
    None
  }
}
