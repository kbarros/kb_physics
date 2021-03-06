package kip.util;

// charAt -> apply

object Parsing {
  def isWhitespace(c: Char): Boolean = {
    c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f'
  }
  
  def stringSplit(ret: Array[String], line: String): Int = {
    var retCnt = 0
    var i = 0
    
    while (i < line.size) {
      while (i < line.size && isWhitespace(line(i))) {
        i += 1
      }
      
      val startIdx = i
      while (i < line.size && !isWhitespace(line(i))) {
        i += 1
      }
      
      if (i > startIdx) {
        ret(retCnt) = line.substring(startIdx, i)
        retCnt += 1
      }
    }
    
    retCnt;
  }
  
  
  def stringSplitDouble(ret: Array[Double], line: String): Int = {
    var retCnt = 0
    var i = 0
    
    while (i < line.size) {
      while (i < line.size && isWhitespace(line(i))) {
        i += 1
      }
      
      val startIdx = i
      while (i < line.size && !isWhitespace(line(i))) {
        i += 1
      }
      
      if (i > startIdx) {
        val str = line.substring(startIdx, i)
        try {
          ret(retCnt) = str.toDouble
          retCnt += 1
        } catch {
          case (nfe: NumberFormatException) => {
            println("Orig line:'"+line+"'")
            println("Indices "+startIdx+","+i)
            println("Failed on str '" + str+"'")
          }
        }
      }
    }
    
    retCnt
  }
}
