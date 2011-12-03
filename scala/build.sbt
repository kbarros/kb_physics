// project name
name := "KB"

// project version
version := "1.0"

// scala version to use for building src
scalaVersion := "2.9.1"

// search subdirectories recursively for unmanaged libraries
unmanagedJars in Compile <++= unmanagedBase map { ub =>
  (ub ** "*.jar").classpath
}

resolvers ++= Seq(
  "Scala-Tools Maven2 Snapshots Repository" at "http://scala-tools.org/repo-snapshots",
  "download.java.net" at "http://download.java.net/maven/2"
)

libraryDependencies ++= Seq(
  "com.twitter" % "json" % "2.1.3",
  "net.liftweb" %% "lift-json" % "2.4-SNAPSHOT",
  "net.java.dev.jna" % "jna" % "3.3.0",
  "org.scala-lang" % "scala-compiler" % "2.9.1",
  "org.scala-lang" % "jline" % "2.9.1"
)

scalacOptions ++= Seq(
  "-deprecation", "-unchecked"
)

javacOptions ++= Seq(
  "-Xlint:unchecked"
)

// Plugin to generate jar file  
seq(sbtassembly.Plugin.assemblySettings: _*)

retrieveManaged := true

// define the statements initially evaluated when entering 'console', 'console-quick', or 'console-project'
initialCommands := """
//  import kip.math.linalg._
//  import DenseMatrix._
  import System.{currentTimeMillis => now}
  def time[T](f: => T): T = {
    val start = now
    try { f } finally { println("Elapsed: " + (now - start)/1000.0 + " s") }
  }
"""

// fork a new JVM for 'run' and 'test:run'
fork := true
