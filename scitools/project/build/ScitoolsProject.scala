import sbt._

class ScitoolsProject(info: ProjectInfo) extends DefaultProject(info)
{
    lazy val hi = task { println("Hello World2  "); None }
    
    def scikit_classpath = "target"/"resources"/"external_classes"
    override def unmanagedClasspath = super.unmanagedClasspath +++ scikit_classpath
}


