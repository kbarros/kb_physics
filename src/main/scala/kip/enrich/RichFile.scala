package kip.enrich
import java.io.File

class RichFile(f: File) {
  def slurp: String = kip.util.Util.readStringFromFile(f.toString)
}
