package kip.util

import kip.util.Util.time
import scala.tools.nsc.interpreter._

object Interpreter {
/*
  def start(bindings: (String, Any)*) {
    val output = new java.io.PrintWriter(new java.io.OutputStreamWriter(Console.out))
    val repl = new scala.tools.nsc.interpreter.ILoop(None, output)
    val settings = new scala.tools.nsc.Settings
    settings.usejavacp.value = true
    settings.classpath.value = System.getProperty("scala.class.path")
    repl.settings = settings
    repl.createInterpreter()
    val varStr = bindings.unzip._1.mkString("[",",","]")
    time(bindings.foreach{ case (k,v) => repl.injectOne(k, v) }, "Binding values "+varStr)
    repl.in = scala.tools.nsc.interpreter.InteractiveReader.createDefault(repl.interpreter)
    try {
      // it is broken on startup; go ahead and exit
      if (repl.interpreter.reporter.hasErrors) return
      repl.printWelcome()
      repl.repl()
    } finally {
      repl.closeInterpreter()
    }
  }
*/
    // copied from scala.tools.nsc.interpreter.ILoop.break

  private def echo(msg: String) = Console println msg
  
  def break[T: Manifest](args: List[NamedParam]): Unit = {
    val msg = if (args.isEmpty) "" else "  Binding " + args.size + " value%s.".format(
      if (args.size == 1) "" else "s"
    )
    echo("Debug repl starting." + msg)
    val repl = new ILoop {
      override def prompt = "\ndebug> "
    }
    repl.settings = new scala.tools.nsc.Settings(echo)
    repl.settings.embeddedDefaults[T]
    repl.createInterpreter()
    repl.in = new JLineReader(new JLineCompletion(repl))

    // rebind exit so people don't accidentally call sys.exit by way of predef
    repl.quietRun("""def exit = println("Type :quit to resume program execution.")""")
    args foreach (p => repl.bind(p.name, p.tpe, p.value))
    repl.loop()

    echo("\nDebug repl exiting.")
    repl.closeInterpreter()
  }      
    

  def main(args: Array[String]) {
//    break(List("a", 1))
    break(Nil)
  }
}
