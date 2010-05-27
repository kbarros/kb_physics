package kip.util;

public class Parsing {

  public static boolean isWhitespace(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f';
  }
  
  
  public static int stringSplit(String[] ret, String line)  {
    int retCnt = 0;
    
    int i = 0;
    while (i < line.length()) {
      while (i < line.length() && isWhitespace(line.charAt(i))) {
        i++;
      }
      
      int startIdx = i;
      while (i < line.length() && !isWhitespace(line.charAt(i))) {
        i++;
      }
      
      if (i > startIdx)
        ret[retCnt++] = line.substring(startIdx, i);
    }
    
    return retCnt;
  }

  
  public static int stringSplitDouble(double[] ret, String line)  {
    int retCnt = 0;
    
    int i = 0;
    while (i < line.length()) {
      while (i < line.length() && isWhitespace(line.charAt(i))) {
        i++;
      }
      
      int startIdx = i;
      while (i < line.length() && !isWhitespace(line.charAt(i))) {
        i++;
      }
      
      if (i > startIdx) {
        String str = line.substring(startIdx, i);
        try {
          double d = Double.parseDouble(str);
          ret[retCnt++] = d; 
        } catch (NumberFormatException nfe) {
          System.out.println("Orig line:'"+line+"'");
          System.out.println("Indices "+startIdx+","+i);
          System.out.println("Failed on str '" + str+"'");
        }
      }
    }
    
    return retCnt;
  }
}
