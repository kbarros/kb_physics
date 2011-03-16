package kip.util

import kip.util.Util.time

object Interpreter {
  def start(bindings: (String, Any)*) {
    val output = new java.io.PrintWriter(new java.io.OutputStreamWriter(Console.out))
    val repl = new scala.tools.nsc.InterpreterLoop(None, output)
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
}
