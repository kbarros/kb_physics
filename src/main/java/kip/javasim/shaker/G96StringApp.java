package kip.javasim.shaker;

import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.sqr;
import static scikit.util.Utilities.format;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.Timer;

import scikit.dataset.Accumulator;
import scikit.dataset.DynamicArray;
import scikit.dataset.Function;
import scikit.dataset.Histogram;
import scikit.graphics.Drawable;
import scikit.graphics.dim2.Gfx2D;
import scikit.graphics.dim2.Scene2D;
import scikit.util.Bounds;
import scikit.util.FileUtil;
import scikit.util.Parsing;
import scikit.util.Terminal;
import scikit.util.Utilities;





public class G96StringApp extends Terminal {
	
	class FrameAnimator {
		FrameDisplay seq;
		Terminal term;
		int delay = 50; // milliseconds
		Timer timer;
		int timerFrame;
		boolean animateStrings = false;
		boolean animateMobile = false;
		int tstar;
		double rstar;
		double rbond;
		File imageDir = null;
		
		
		public FrameAnimator(FrameDisplay seq) {
			this.seq = seq;
			this.term = seq.term;
			timer = new Timer(delay, taskPerformer);
		}
		
		ActionListener taskPerformer = new ActionListener() {
		    public void actionPerformed(ActionEvent evt) {
		    	if (timerFrame >= seq.nframes) {
		    		timer.stop();
		    	}
		    	else {
		    		if (animateMobile)
		    			seq.drawMobile(timerFrame, tstar, rstar);
		    		else if (animateStrings)
		    			seq.drawStrings(timerFrame, tstar, rstar, rbond);
		    		else
		    			seq.draw(timerFrame);
		    		if (imageDir != null) {
						try {
							String fname = Utilities.formatI4(timerFrame);
							String imgName = imageDir+File.separator+fname+".png";
							ImageIO.write(seq.plot.getImage(320, 320), "png", new File(imgName));
						} catch (IOException e) {
							term.println(e.getMessage());
						}
		    		}
		    		timerFrame++;
		    	}
		    }
		};
		
		void setMobileParams(int tstar, double rstar) {
			animateMobile = true;
			this.tstar = tstar;
			this.rstar = rstar;
			timerFrame = max(tstar+2, timerFrame);
		}

		void setStringParams(int tstar, double rstar, double rbond) {
			animateStrings = true;
			this.tstar = tstar;
			this.rstar = rstar;
			this.rbond = rbond;
			timerFrame = max(tstar+2, timerFrame);
		}

		public void start() {
			timer.start();
		}
		
		public void stop() {
			if (timer != null)
				timer.stop();
		}
		
		public void selectImageDir() {
			try {
				imageDir = FileUtil.directoryDialog("Select output directory");
				term.println("Directory: " + imageDir);
			} catch(IOException e) {
				term.println(e);
			}
		}
	}

	class FrameDisplay {
		Terminal term;
		ParticleData data;
		double radius;
		int nframes;
		Bounds bds;
		Scene2D plot = new Scene2D("Atoms");
		
		public FrameDisplay(Terminal term, ParticleData data, double radius) {
			this.term = term;
			this.data = data;
			this.radius = radius;
			nframes = data.frames.size();
			bds = new Bounds(0, data.xw, 0, data.yw, 0, data.zw);
			draw(0);
			Utilities.frame(plot.getComponent(), plot.getTitle());
		}
		
		public FrameAnimator animator() {
			return new FrameAnimator(this);
		}
		
		
		public Drawable<Gfx2D> getBackgroundDrawable(final int t) { 
			return new Drawable<Gfx2D>() {
				public void draw(Gfx2D g) {
					g.setColor(Color.BLACK);
					g.drawRect(bds.xmin, bds.ymin, bds.xmax-bds.xmin, bds.ymax-bds.ymin);
					
					g.setProjection(g.pixelBounds());
					
					g.setColor(new Color(1f, 1f, 1f, 0.5f));
					g.fillRect(0, 0, g.pixelBounds().getWidth(), 22);
					
					g.setColor(Color.RED);
					int particles = data.nparticles;
					String s = "Frame="+t+", NumAtoms="+particles+", BoxWidth=<"+format(data.xw)+","+format(data.yw)+","+format(data.zw)+">";
					g.drawString(s, 4, 4);
					
					g.setProjection(g.viewBounds());
				}
				public Bounds getBounds() {
					return bds;
				}
			};
		}
		

		public Scene2D draw(int t) {
			if (t >= nframes)
				term.println("Time " + t + " exceeds maximum frame " + (nframes-1));
			plot.clearDrawables();
			plot.addDrawable(data.frames.get(t).getDrawable(new Color(0f, 0f, 1f, 0.5f), radius));
			plot.addDrawable(getBackgroundDrawable(t));
			return plot;
		}
		
		public void drawMobile(int t, int tstar, double rstar) {
			if (!data.sanityCheckFrameNumbers(t, tstar))
				return;

			Frame mobile = data.getMobileParticles(t, tstar, rstar);
			plot.clearDrawables();
			plot.addDrawable(data.frames.get(t).getDrawable(new Color(0f, 0f, 1f, 0.2f), radius));
			plot.addDrawable(mobile.getDrawable(new Color(0f, 1f, 0.2f, 1f), radius));
			plot.addDrawable(getBackgroundDrawable(t));
		}
		
		
		public void drawStrings(int t, int tstar, double rstar, double rbond) {
			drawStrings(t, data.getStrings(t, tstar, rstar, rbond));
		}
		
		public void drawStrings(int t, List<Frame> strings) {
			if (t >= nframes)
				term.println("Time " + t + " exceeds maximum frame " + (nframes-1));

			plot.clearDrawables();
			plot.addDrawable(data.frames.get(t).getDrawable(new Color(0f, 0f, 1f, 0.2f), radius));
			
			Color colors[] = {Color.CYAN, Color.GRAY, Color.MAGENTA, Color.ORANGE, Color.PINK, Color.RED};
			for (int i = 0; i < strings.size(); i++) {
				Frame f = strings.get(i);
				plot.addDrawable(f.getDrawable(colors[i%colors.length], radius));
			}
			plot.addDrawable(getBackgroundDrawable(t));
		}
	}

	class Frame {
		int t;
		public DynamicArray id = new DynamicArray();
		public DynamicArray x = new DynamicArray();
		public DynamicArray y = new DynamicArray();
		public DynamicArray z = new DynamicArray();
		
		public Frame(int t) {
			this.t = t;
		}
		
		public void append(double id, double x, double y, double z) {
			this.id.append(id);
			this.x.append(x);
			this.y.append(y);
			this.z.append(z);
		}
		
		public Drawable<Gfx2D> getDrawable(final Color color, final double radius) {
			return new Drawable<Gfx2D>() {
				public void draw(Gfx2D g) {
					g.setColor(color);
					for (int i = 0; i < id.size(); i++) {
						g.fillCircle(x.get(i), y.get(i), radius);
					}
				}
				public Bounds getBounds() {
					return new Bounds();
				}
			};
		}
	}

	class ParticleData {
		Terminal term;
		public int frameJump = 1;
		ArrayList<Frame> frames = new ArrayList<Frame>();
		double xw, yw, zw;
		int nparticles;
		
		
		boolean stringsMatch(String s1, String s2) {
			return s1 != null && s2 != null & s1.trim().equals(s2.trim());
		}
		
		String verifyEntry(String received, String expected) {
			if (received == null) {
				System.err.println("Got null instead of " + expected);
			}
			if (!stringsMatch(received, expected)) {
				System.err.println("Error parsing on line: "+expected+" "+received);
			}
			return received;
		}
		
		
//		void parse1(Terminal term, String fname, int maxFrames) throws FileNotFoundException {
//			Scanner sc = new Scanner(new File(fname));
//			
//			verifyEntry(sc.nextLine(), "TITLE");
//			term.println("File description: "+sc.nextLine());
//			verifyEntry(sc.nextLine(), "END");
//			
//			term.print("Reading frames... ");
//			for (int t = 0; t < maxFrames && sc.hasNext("TIMESTEP"); t++) {
//				// read time step
//				verifyEntry(sc.nextLine(), "TIMESTEP");
//				sc.nextLine(); // throw away
//				verifyEntry(sc.nextLine(), "END");
//				
//				// read particle coordinates
//				verifyEntry(sc.nextLine(), "POSITIONRED");
//				Frame f = new Frame(t);
//				while (sc.hasNextDouble()) {
//					f.x.append(sc.nextDouble());
//					f.y.append(sc.nextDouble());
//					f.z.append(sc.nextDouble());
//					sc.nextLine();
//				}
//				frames.add(f);
//				verifyEntry(sc.nextLine(), "END");
//				if (t == 0) {
//					nparticles = f.x.size();
//				}
//				else {
//					if (nparticles != f.x.size())
//						term.println("Error: inconsistent particle numbers "+nparticles+" and "+f.x.size());
//				}
//				
//				// read box size
//				verifyEntry(sc.nextLine(), "BOX");
//				if (t == 0) {
//					xw = sc.nextDouble();
//					yw = sc.nextDouble();
//					zw = sc.nextDouble();
//					sc.nextLine();
//				}
//				else {
//					if (xw != sc.nextDouble() || yw != sc.nextDouble() || zw != sc.nextDouble())
//						term.println("Error: Box size non-constant.");
//					sc.nextLine();
//				}
//				verifyEntry(sc.nextLine(), "END");
//				
//				if (t % 100 == 0)
//					term.print(t + " ");
//			}
//			
//			term.println("done. "+ frames.size() + " frames read.");
//		}
		
		
		void parse2(String fname, int maxFrames, int skipEvery) throws IOException {
			@SuppressWarnings("resource")
            BufferedReader is = new BufferedReader (new FileReader (fname));
			
			verifyEntry(is.readLine(), "TITLE");
			term.println("File description: "+is.readLine());
			verifyEntry(is.readLine(), "END");
			
			term.print("Reading frames... ");
			
			int t;
			for (t = 0; t < maxFrames && stringsMatch(is.readLine(), "TIMESTEP"); t++) {
				is.readLine(); // throw away time step information
				verifyEntry(is.readLine(), "END");
				
				// read particle coordinates
				verifyEntry(is.readLine(), "POSITIONRED");
				Frame f = new Frame(t);
				int particleCnt = 0;
				double[] pos = new double[3];
				for (String line = is.readLine(); !line.equals("END"); line = is.readLine()) {
					int cnt = Parsing.stringSplitDouble(pos, line);
					if (cnt != 3) {
						System.err.println("Cannot parse position from line: '" + line +"'");
					}
					f.append(particleCnt++, pos[0], pos[1], pos[2]);
				}
				
				if (t % skipEvery == 0)
					frames.add(f);
				
				if (t == 0) {
					nparticles = particleCnt;
				}
				else {
					if (nparticles != particleCnt)
						term.println("Error: inconsistent particle numbers "+nparticles+" and "+particleCnt);
				}
				
				// read box size
				verifyEntry(is.readLine(), "BOX");
				String line = is.readLine();
				int cnt = Parsing.stringSplitDouble(pos, line);
				if (cnt != 3) {
					System.err.println("Cannot parse box size from line: '" + line +"'");
				}
				if (t == 0) {
					xw = pos[0];
					yw = pos[1];
					zw = pos[2];
				}
				else {
					if (xw != pos[0] || yw != pos[1] || zw != pos[2])
						term.println("Error: Box size non-constant.");
				}
				verifyEntry(is.readLine(), "END");
				
				if (t % 500 == 0)
					term.print(t + " ");
			}
			
			term.println("done. Read "+ frames.size() + " of " + t + " frames.");
			term.println("Particle number: " + nparticles);
		}
		
		
		public ParticleData(Terminal term, String fname, int maxFrames, int skipEvery) {
			this.term = term;
			
			try {
				parse2(fname, maxFrames, skipEvery);
			} catch (IOException e) {
				term.println(e);
			}
		}
		
		public Frame getMobileParticles(int t, int tstar, double rstar) {
			if (!sanityCheckFrameNumbers(t, tstar))
				return null;
			
			int t2 = t;
			int t1 = t2 - tstar;

			Frame f1 = frames.get(t1);
			Frame f2 = frames.get(t2);
			Frame mobile = new Frame(t2);

			for (int i = 0; i < nparticles; i++) {
				double x1 = f1.x.get(i);
				double y1 = f1.y.get(i);
				double z1 = f1.z.get(i);
				double x2 = f2.x.get(i);
				double y2 = f2.y.get(i);
				double z2 = f2.z.get(i);

				double d2 = dist2(x2-x1, y2-y1, z2-z1);
				if (d2 > sqr(rstar)) {
					mobile.append(i, x2, y2, z2);
				}
			}
			return mobile;
		}
		
		public boolean sanityCheckFrameNumbers(double t, double tstar) {
			if (t >= frames.size()) {
				term.println("ERROR: Frame " + t + " exceeds maximum, " + (frames.size()-1));
				return false;
			}
			if (t-tstar < 0) {
				term.println("ERROR: Frame " + t + " minus tstar="+tstar+" is too early");
				return false;
			}
			return true;
		}
		
		double dist2(double dx, double dy, double dz) {
			while (dx >  xw/2.) dx -= xw;
			while (dx < -xw/2.) dx += xw;
			while (dy >  yw/2.) dy -= yw;
			while (dy < -yw/2.) dy += yw;
			while (dz >  zw/2.) dz -= zw;
			while (dz < -zw/2.) dz += zw;
			return dx*dx + dy*dy + dz*dz;
		}
		
		
		public List<Frame> getStrings(int t, int tstar, double rstar, double rbond) {
			if (!sanityCheckFrameNumbers(t, tstar))
				return null;

			Frame mobile = getMobileParticles(t, tstar, rstar);
			
			Percolator<Integer> perc = new Percolator<Integer>();
			for (int i = 0; i < mobile.x.size(); i++) {
				Frame f_pre = frames.get(t-tstar); 
				if (f_pre == null)
					System.err.println("failed on frame " + (t-tstar));
				
				double x1_cur = mobile.x.get(i);
				double y1_cur = mobile.y.get(i);
				double z1_cur = mobile.z.get(i);
				double x1_pre = f_pre.x.get(i);
				double y1_pre = f_pre.y.get(i);
				double z1_pre = f_pre.z.get(i);
				
				perc.add(i);
				
				for (int j = 0; j < i; j++) {
					double x2_cur = mobile.x.get(j);
					double y2_cur = mobile.y.get(j);
					double z2_cur = mobile.z.get(j);
					double x2_pre = f_pre.x.get(j);
					double y2_pre = f_pre.y.get(j);
					double z2_pre = f_pre.z.get(j);
					
					if (dist2(x1_cur-x2_pre, y1_cur-y2_pre, z1_cur-z2_pre) < sqr(rbond) ||
						dist2(x1_pre-x2_cur, y1_pre-y2_cur, z1_pre-z2_cur) < sqr(rbond))
						perc.bond(i, j);
				}
			}
			
			ArrayList<Frame> ret = new ArrayList<Frame>();
			for (ArrayList<Integer> group : perc.getGroups()) {
				Frame f = new Frame(t);
				for (int i : group) {
					f.append(i, mobile.x.get(i), mobile.y.get(i), mobile.z.get(i));
				}
				ret.add(f);
			}
			return ret;
		}
	}

	class Alpha {
		public Accumulator r2;
		public Accumulator r4;
		public Accumulator alpha;
	}

	
	class Commands {
		Terminal term;
		
		public Commands(Terminal term) {
			this.term = term;
		}
		
//		public ParticleData loadData(String fname) {
//			term.println("Opening " + fname);
//			return new ParticleData(term, fname);
//		}
		
		public ParticleData loadData() {
			return loadData(Integer.MAX_VALUE, 1);
		}
		
		public ParticleData loadData(int maxFrames, int skipEvery) {
			try {
				String fname = FileUtil.loadDialog(term.getConsole(), "Open G96 Particle Data");
				if (fname == null)
					return null;
				else {
					term.println("Opening " + fname);
					return new ParticleData(term, fname, maxFrames, skipEvery);
				}
			} catch(IOException e) {
				term.println(e);
				return null;
			}
		}
		
		public Alpha diffusionAnalysis(ParticleData data) {
			int nsteps = 500;
			int[] steps = new int[nsteps];
			for (int i = 0; i < nsteps; i++) {
				steps[i] = (int)exp(i*0.1);
			}
			
			double bw = 1;
			Accumulator r2 = new Accumulator(bw);
			Accumulator r4 = new Accumulator(bw);

			for (int t1 = 0; t1 < data.frames.size(); t1 += data.frameJump) {
				Frame f1 = data.frames.get(t1);
				for (int s = 0; s < nsteps; s++) {
					int t2 = t1+steps[s];
					if (t2 >= data.frames.size()) break;
					Frame f2 = data.frames.get(t2);

					for (int i = 0; i < data.nparticles; i++) {
						double dx = f2.x.get(i)-f1.x.get(i);
						double dy = f2.y.get(i)-f1.y.get(i);
						double dz = f2.z.get(i)-f1.z.get(i);
						double d2 = data.dist2(dx, dy, dz);
						r2.accum(t2-t1, d2);
						r4.accum(t2-t1, d2*d2);
					}
				}
			}
			
			Accumulator alpha = new Accumulator(bw);
			for (double t : r2.keys()) {
				double dim = 3;
				double d2 = r2.eval(t);
				double d4 = r4.eval(t);
				alpha.accum(t, (1/(1.+2./dim)) * (d4/(d2*d2)) - 1.0);
			}

			Alpha ret = new Alpha();
			ret.r2 = r2;
			ret.r4 = r4;
			ret.alpha = alpha;
			return ret;
		}
		
		public Histogram vanHove(ParticleData data, int tstar) {
			Histogram ret = new Histogram(0.05);
			ret.setNormalizing(true);
			
			for (int t1 = 0; t1 < data.frames.size()-tstar; t1 += data.frameJump) {
				int t2 = t1 + tstar;
				Frame f1 = data.frames.get(t1);
				Frame f2 = data.frames.get(t2);

				for (int i = 0; i < data.nparticles; i++) {
					double dx = f2.x.get(i)-f1.x.get(i);
					double dy = f2.y.get(i)-f1.y.get(i);
					double dz = f2.z.get(i)-f1.z.get(i);
					double r = sqrt(data.dist2(dx, dy, dz));
					ret.accum(r, r);
				}
			}
			return ret;
		}
		
		public Function vanHoveTheory(final double sig2) {
			return new Function(0, 5*sqrt(sig2)) {
				public double eval(double r) {
					return (4*PI*r*r) * Math.pow(1/(sig2*2*PI), 3./2) * exp(-(r*r)/(2*sig2));
				}
			};
		}
		
		public double countMobile(ParticleData data, int tstar, double rstar) {
			double mobile = 0;
			int nframes = 0;
			for (int t1 = tstar+2; t1 < data.frames.size()-tstar; t1 += data.frameJump) {
				mobile += data.getMobileParticles(t1, tstar, rstar).x.size();
				nframes++;
			}
			return mobile / nframes;
		}
		
		public FrameDisplay makeDisplay(ParticleData data, double radius) {
			return new FrameDisplay(term, data, radius);
		}
	}

	
	public static void main(String[] args) {
		G96StringApp stringApp = new G96StringApp();
		stringApp.help = "Suggested commands:\n"+
			"(Plots are dynamic: Right click to change view options or to save data,\n"+
			" left click and drag to zoom in, double click to zoom out.)\n"+
			"\n"+
			"\tdata = loadData(100000, 10); // analyze only 1 in 10 of the first 100,000 frames\n"+
			"\tdata = loadData(); // or analyze all the data (requires a lot of RAM)\n"+
			"\n"+
			"\t// if the particles diffuse slowly, we can speed diffusionAnalysis() and\n" +
//			"\t// vanHove() by skipping 'frameJump-1' frames for each one to be analyzed\n"+
//			"\tdata.frameJump = 100; // default 100\n\n"+
			"\ta = diffusionAnalysis(data);\n"+
			"\tplot(a.r2, \"<r^2> vs frame#\"); // mean squared displacement\n"+
			"\tplot(a.alpha, \"alpha vs frame#\"); // non-Gaussian parameter 3<r^4>/5<r^2>^2\n"+
			"\t// (find t* from maximum of alpha)\n\n"+
			"\t// plot van-hove self correlation using, e.g., t* = 300\n"+
			"\tplot(vanHove(data, 300), \"Data\");\n"+
			"\t// compare to Gaussian using, e.g., sigma^2=0.05 (the value measured from <r^2(t*)>)\n"+
			"\treplot(vanHoveTheory(0.05), \"Theory\");\n"+
			"\t// (find r* from the right-most intersection of the two plots)\n"+
			"\n"+
			"\t// Display the data\n"+
			"\ts = makeDisplay(data, 0.03); // visualize particles using circles of radius 0.03\n"+
			"\ts.draw(100); // draw the 101th frame of the data\n"+
			"\t// draw \"mobile\" particles using at frame 1000 using, e.g., t*=300, r*=0.7\n"+
			"\ts.drawMobile(1000, 300, 0.7);\n"+
			"\n"+
			"\t// String analysis\n"+
			"\t// Get a list of strings using at frame 1000 using, e.g., t*=300, r*=0.7, rbond=0.1\n"+
			"\t// (Choose rbond to be less than the hard core diameter)\n"+
			"\tstrs = data.getStrings(1000, 300, 0.7, 0.1);\n"+
			"\ts.drawStrings(1000, strs); // plot strings at frame 1000\n"+
			"\n"+
			"\t// Animate the data\n"+
			"\ta = s.animator();\n" +
			"\ta.start(); // animate the tracked particles\n"+
			"\ta.stop();\n" +
			"\ta.setMobileParams(300, 0.7); // set (t*, r*) for mobile particles\n"+
			"\ta.selectImageDir(); // (optional) write animation frames to images\n"+
			"\ta.start(); // animate mobile particles\n"+
			"";
		stringApp.importObject(stringApp.new Commands(stringApp));
		stringApp.runApplication();
	}
}
