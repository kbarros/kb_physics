package kip.javasim.shaker;

import java.io.*;
//import java.util.Scanner;
//import java.nio.ByteBuffer;
//import java.nio.MappedByteBuffer;
//import java.nio.channels.FileChannel;


public class ScanTestApp {
	
	static boolean isWhitespace(char c) {
		return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f';
	}
	
	
	static int fastStringSplit(String[] ret, String line)  {
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
			ret[retCnt++] = line.substring(startIdx, i);
		}
		
		return retCnt;
	}

	
	static int fastStringSplitDouble(double[] ret, String line)  {
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
			
			String str = line.substring(startIdx, i);
			
            try {
                double d = Double.parseDouble(str);
    			ret[retCnt++] = d; 
            } catch (NumberFormatException nfe) {
            	System.out.println("Failed on str " + str);
            }
		}
		
		return retCnt;
	}

	static void streamFile2(File file) {
		try {
			@SuppressWarnings("resource")
			BufferedReader is = new BufferedReader (new FileReader (file));
//			String[] splitStrs = new String[256];
			double[] splitVals = new double[256];
			
			long i = 0;
			while (true) {
				String b = is.readLine();
				if (b == null)
					break;
				i += b.length();
				int cnt = fastStringSplitDouble(splitVals, b);
				i += cnt / 2000;
//				fastStringSplitDouble(splitVals, b);
//				i += splitVals.length / 2000;
				
				long Mb = i / (1024L*1024L);
				if (i % (10L*1024L*1024L) == 0)
					System.out.println(Mb + "MB "  +b);
				
//				if (Mb >= 200) break;
			}
			long Mb = i / (1024L*1024L);
			System.out.println(Mb + "MB ");
		}
		catch (IOException e) {
			System.out.println("file not found "  + e);
		}

	}

	
	
//	static void streamFile(File file) {
//		try {
//			InputStream is = new BufferedInputStream (new FileInputStream (file));
//			
//			for (long i = 0; i < Integer.MAX_VALUE; i++) {
//				int b = is.read();
//				if (b == -1)
//					break;
//				
//				long Mb = i / (1024L*1024L);
//				if (i % (10L*1024L*1024L) == 0)
//					System.out.println(Mb + "MB "  +b);
//				
//				if (Mb >= 200) break;
//			}
//		}
//		catch (IOException e) {
//			System.out.println("file not found "  + e);
//		}
//
//	}
//
//
//	static void scanFile(File file) {
//		try {
//			Scanner sc = new Scanner(file);
//			long bytes = 0;
//			while (sc.hasNextLine()) {
//				String line = sc.nextLine();
//				bytes += line.length();
//				
//				long Mb = bytes / (1024L*1024L);
//				if (bytes % (10L*1024L*1024L) == 0) {
//					System.out.println(Mb + "MB " + line);
//				}
//				
//				if (Mb >= 200) break;
//			}
//		}
//		catch (FileNotFoundException e) {
//			System.out.println("file not found "  + e);
//		}
//
//	}
//	
//	static void scanFile2(File file) {
//		try {
//			Scanner sc = new Scanner(file);
//			long i = 0;
//			while (sc.hasNextLine()) {
////				if (sc.hasNextDouble())
////					sc.nextDouble();
////				else {
//				sc.nextLine();
//				for (int j = 0; j < 10; j++)
//					Double.parseDouble(Double.toString(j));
//				
//				i++;
//				if (i % 100000 == 0)
//					System.out.println(i);
//			}
//		}
//		catch (FileNotFoundException e) {
//			System.out.println("file not found "  + e);
//		}
//
//	}

//	private static int sum(ByteBuffer bb) {
//		int sum = 0;
//		while (bb.hasRemaining()) {
//			if ((sum & 1) != 0)
//				sum = (sum >> 1) + 0x8000;
//			else
//				sum >>= 1;
//				sum += bb.get() & 0xff;
//				sum &= 0xffff;
//		}
//		return sum;
//	}
//	static void mapFile(File file) {
//	    try {
//	    	// Open the file and then get a channel from the stream
//	    	FileInputStream fis = new FileInputStream(file);
//	    	FileChannel fc = fis.getChannel();
//
//	    	// Get the file's size and then map it into memory
//	    	long sz = fc.size();
//	    	System.out.println(sz);
//	    	sz = Integer.MAX_VALUE;
//	    	MappedByteBuffer bb = fc.map(FileChannel.MapMode.READ_ONLY, 0, sz);
//
//	    	// Compute and print the checksum
//	    	int sum = sum(bb);
//	    	long M =1024L*1024L; 
//	    	long Mb = (sz + (M-1)) / M;
//	    	String s = Integer.toString(sum);
//	    	System.out.println(s + "\t" + Mb + "\t" + file);
//
//	    	// Close the channel and the stream
//	    	fc.close();
//	    } catch (IOException e) {
//			System.err.println(file + ": " + e);
//	    }
//	}
	
	public static void main(String[] args) {
		System.out.println(System.getProperty("java.version"));

		File file = new File("/Users/kbarros/Desktop/data/HIVmono_10000ps.g96");
		long t1, t2;
		
//		t1 = System.currentTimeMillis();
//		System.out.println("begin scan");
//		scanFile2(file);
//		t2 = System.currentTimeMillis();
//		System.out.println("end scan " + ((t2-t1))/1000. + "s");
		
//		t1 = System.currentTimeMillis();
//		System.out.println("begin map");
//		mapFile(file);
//		t2 = System.currentTimeMillis();
//		System.out.println("end map " + ((t2-t1))/1000. + "s");

		t1 = System.currentTimeMillis();
		System.out.println("begin stream");
		streamFile2(file);
		t2 = System.currentTimeMillis();
		System.out.println("end stream " + ((t2-t1))/1000. + "s");

	}
	
}
