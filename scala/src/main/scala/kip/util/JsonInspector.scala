package kip.util

case class JsonInspector(tree: Any) {
  def toDouble: Double = {
    tree match {
      case v:Int    => v.toDouble
      case v:Double => v
      case v:BigDecimal => v.toDouble
      case v:AnyRef => sys.error("Cannot convert class "+v.getClass+" to double")
    }
  }
  
  def toInt: Int = {
    tree match {
      case v:Int => v
    }
  }
  
  def apply(key: Any): JsonInspector = {
    tree match {
      case m: Map[Any,Any] => JsonInspector(m(key))
    }
  }
  
  def isEmpty: Boolean = {
    tree match {
      case m: Map[Any, Any] => m.isEmpty
      case l: List[Any] => l.isEmpty
      case _ => false
    }
  }
  
  def ++(that: JsonInspector): JsonInspector = {
    (tree, that.tree) match {
      case (m1: Map[Any, Any], m2: Map[Any, Any]) => JsonInspector(m1++m2)
    }
  }
  
  
}
