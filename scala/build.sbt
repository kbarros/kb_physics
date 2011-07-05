// project name
name := "KB"

// project version
version := "1.0"

// scala version to use for building src
scalaVersion := "2.9.0"

// scala version for building project definition ??


libraryDependencies ++= Seq(
  "net.java.dev.jna" % "jna" % "3.3.0",
  "com.twitter" % "json" % "2.1.3",
  "org.scala-lang" % "scala-compiler" % "2.9.0",
  "org.scala-lang" % "jline" % "2.9.0"
)

resolvers ++= Seq(
  "Scala-Tools Maven2 Snapshots Repository" at "http://scala-tools.org/repo-snapshots",
  "download.java.net" at "http://download.java.net/maven/2"
)

scalacOptions ++= Seq(
  "-deprecation", "-unchecked"
)

javacOptions ++= Seq(
  "-Xlint:unchecked"
)

retrieveManaged := true

// define the statements initially evaluated when entering 'console', 'console-quick', or 'console-project'
initialCommands := """
  import kip.math.linalg._
  import DenseMatrix._
/*
  import scalala.scalar._
  import scalala.tensor.::
  import scalala.tensor.mutable._
  import scalala.tensor.dense._
  import scalala.tensor.sparse._
  import scalala.library.Library._
  import scalala.library.LinearAlgebra._
  import scalala.library.Statistics._
  import scalala.library.Plotting._
  import scalala.operators.Implicits._
*/
  import System.{currentTimeMillis => now}
  def time[T](f: => T): T = {
    val start = now
    try { f } finally { println("Elapsed: " + (now - start)/1000.0 + " s") }
  }
"""

// fork a new JVM for 'run' and 'test:run'
fork := true
