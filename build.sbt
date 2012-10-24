name := "KB"

version := "1.0"

scalaVersion := "2.9.1"

// search subdirectories recursively for unmanaged libraries
unmanagedJars in Compile <++= unmanagedBase map { ub =>
  (ub ** "*.jar").classpath
}

resolvers ++= Seq(
  "download.java.net" at "http://download.java.net/maven/2"
)

libraryDependencies ++= Seq(
  "com.twitter" % "json" % "2.1.3",
  "net.liftweb" %% "lift-json" % "2.4",
  "net.java.dev.jna" % "jna" % "3.3.0",
  "org.scala-lang" % "scala-compiler" % "2.9.1",
  "org.scala-lang" % "jline" % "2.9.1",
  "com.googlecode.matrix-toolkits-java" % "mtj" % "0.9.14"
)

scalacOptions ++= Seq(
  "-deprecation", "-unchecked"
)

javacOptions ++= Seq(
  "-Xlint:unchecked"
)

retrieveManaged := true

// fork a new JVM for 'run' and 'test:run'
fork := true
