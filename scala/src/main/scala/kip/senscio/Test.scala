package kip.senscio

/*
object Test {
  
  def testLoad() {
    val file = scala.io.Source.fromFile("/Users/kbarros/Desktop/FinalTest.txt")
    val lines = file.getLines
    val first = lines.next()
    val keys = first.split("\t") 
    
    for (line <- lines) {
      val data = line.split("\t")
      val withdrawal = data(3)
      val name = data(5)
      println(withdrawal + " " + name)
      
    }
  }
  
  
}



object Roll {
  
  def simpleImpl: Roll = new Roll {
    
    var classes = Map[String, Clazz]()
    var concepts = Seq[Concept]()
    var links = Map[Concept, Seq[Link]]()
    
    
    def defineClass(name: String, parent: Option[Clazz]): Clazz = {
      val ret = new Clazz {
        def name = name
        def parent = parent
      }
      classes += (name -> ret)
      ret
    }
    
    def classForName(name: String) {
      classes(name)
    }
  
    def createConcept(clazz: Clazz, data: Any): Concept = {
      val ret = new Concept {
        def ofClass = clazz
        def data = data
      }
      ret
    }
    
    def createLink(subj: Concept, obj: Concept) {
      links.getOrElse(subj, 
    }
  
  def processDataFile(lines: Traversable[String], ctors: Seq[String => Clazz])
  
  def linksForConcept(concept: Concept): Seq[Concept]
  
  
  def triggerOnConceptCreation()

  }
}

abstract class Roll {
  
  
  trait Clazz {
    def name: String
    def parent: Option[Clazz]
  }
  
  trait Concept {
    def ofClass: Clazz
    def data: Any
  }
  
  trait Link {
    def getSubject: Concept
    def getObject: Concept
  }
  
  def defineClass(name: Symbol, parent: Option[Clazz]): Clazz
  def classForName(name: Symbol)
  
  def createConcept(clazz: Clazz, data: Any): Concept
  def createLink(subj: Concept, obj: Concept)
  
  def processDataFile(lines: Traversable[String], ctors: Seq[String => Clazz])
  
  def linksForConcept(concept: Concept): Seq[Concept]
  
  
  def triggerOnConceptCreation()
}
*/
