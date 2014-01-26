name := "KB"

version := "1.0"

scalaVersion := "2.10.3"

// search subdirectories recursively for unmanaged libraries
unmanagedJars in Compile <++= unmanagedBase map { ub =>
  (ub ** "*.jar").classpath
}

resolvers ++= Seq(
  "download.java.net" at "http://download.java.net/maven/2"
)

libraryDependencies ++= Seq(
  "com.fasterxml.jackson.module" %% "jackson-module-scala" % "2.3.0",
  "net.java.dev.jna" % "jna" % "3.3.0",
  "org.scala-lang" % "scala-compiler" % "2.10.3",
  "org.scala-lang" % "jline" % "2.10.3",
  "com.googlecode.matrix-toolkits-java" % "mtj" % "1.0.1" // "0.9.14"
  // "com.github.fommil.netlib" % "all" % 1.1.2
)

scalacOptions ++= Seq(
  "-deprecation", "-unchecked",
  "-language:implicitConversions", "-language:higherKinds"
)

javacOptions ++= Seq(
  "-Xlint:unchecked"
)

retrieveManaged := true

// fork a new JVM for 'run' and 'test:run'
fork := true
