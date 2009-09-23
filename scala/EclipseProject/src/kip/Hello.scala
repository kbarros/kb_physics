package kip

import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation


object Hello {
	def main(args : Array[String]) {
		val x = "yo"
		println ("hello" + x)
		ActorTest.go()
		new Control(new Hello(), "Hello")
	}
}


class Hello extends Simulation {
	def clear() {
		println ("clear")
	}
	def animate() {
		println ("animating")
	}
	def run() {
		var a = List()
		while (true) {
			println ("running")
			Job.animate()
		}
	}
	
	override def load(c: Control) {
//		c.frame(histogramPlot, densityPlot, heatPlot);
//		params.add("L", 16);
//		params.add("Impurity density", 0.2);
//		params.add("mcs");
//		params.add("iteration");
	}

}

import scala.actors.Actor
import scala.actors.Actor._


object ActorTest {
	object Ping
	object Pong
	object Stop


	def pinger(count: Int, pong: Actor) = actor {
			var pingsLeft = count - 1
			pong ! Ping
			loop {
				react {
				case Pong =>
				if (pingsLeft % 10000 == 0)
					Console.println("Ping: pong")
					if (pingsLeft > 0) {
						pong ! Ping
						pingsLeft -= 1
					} else {
						Console.println("Ping: stop")
						pong ! Stop
						exit()
					}
				}
			}
	}
	
	val pong = actor  {
		var pongCount = 0
		loop {
			react {
			case Ping =>
			if (pongCount % 10000 == 0)
				Console.println("Pong: ping "+pongCount)
				sender ! Pong
				pongCount = pongCount + 1
			case Stop =>
			Console.println("Pong: stop")
			exit()
			}
		}
	}


	def go() {
		val ping = pinger(1000000, pong)
		ping.start
		pong.start
	}
}
