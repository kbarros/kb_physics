package kip.util

import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.module.scala.DefaultScalaModule
import com.fasterxml.jackson.module.scala.experimental.ScalaObjectMapper

object JacksonWrapper {
  val mapper = new ObjectMapper() with ScalaObjectMapper 
  mapper.registerModule(DefaultScalaModule)

  def serialize(value: Any): String = {
    mapper.writeValueAsString(value)
  }
  
  def deserialize[T: Manifest](value: String) : T =
    mapper.readValue[T](value)
}
