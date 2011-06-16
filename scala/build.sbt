// project name
name := "KB"

// project version
version := "1.0"

// scala version to use for building src
scalaVersion := "2.9.0"

// scala version for building project definition ??


libraryDependencies ++= Seq(
  "net.java.dev.jna" % "jna" % "3.2.3",
  "com.twitter" % "json" % "2.1.3",
  "org.scalala" %% "scalala" % "1.0.0.RC2-SNAPSHOT",
  "org.scala-lang" % "scala-compiler" % "2.9.0"
)

resolvers ++= Seq(
  "Scala-Tools Maven2 Snapshots Repository" at "http://scala-tools.org/repo-snapshots",
  "ScalaNLP" at "http://repo.scalanlp.org/repo/",
  "Ondex" at "http://ondex.rothamsted.bbsrc.ac.uk/nexus/content/groups/public/"
)

scalacOptions ++= Seq(
  "-deprecation", "-unchecked"
)

retrieveManaged := true


// define the statements initially evaluated when entering 'console', 'console-quick', or 'console-project'
initialCommands := """
  import System.{currentTimeMillis => now}
  def time[T](f: => T): T = {
    val start = now
    try { f } finally { println("Elapsed: " + (now - start)/1000.0 + " s") }
  }
"""

// fork a new JVM for 'run' and 'test:run'
fork := true
