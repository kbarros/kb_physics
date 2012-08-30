package kip.javasim.shaker;

import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.sqr;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Hashtable;
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
import scikit.util.Terminal;
import scikit.util.Utilities;

class FrameAnimator {
	FrameSequence seq;
	Terminal term;
	int delay = 50; // milliseconds
	Timer timer;
	int timerFrame;
	boolean animateMobile = false;
	int tstar;
	double rstar;
	File imageDir = null;
	
	
	public FrameAnimator(FrameSequence seq) {
		this.seq = seq;
		this.term = seq.term;
		timer = new Timer(delay, taskPerformer);
	}
	
	ActionListener taskPerformer = new ActionListener() {
	    public void actionPerformed(ActionEvent evt) {
	    	if (timerFrame >= seq.frames.length) {
	    		timer.stop();
	    	}
	    	else {
	    		if (!animateMobile)
	    			seq.plot(timerFrame);
	    		else
	    			seq.plotMobile(timerFrame, tstar, rstar);
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

class FrameSequence {
	Terminal term;
	Frame[] frames;
	int nframes;
	Bounds bds;
	Scene2D plot = new Scene2D("Tracked particles");
	
	public FrameSequence(Terminal term, ShakerData d) {
		this.term = term;
		
		term.print("Counting frames... ");
		calculateSceneBounds(d);
		term.println(nframes + " counted.");
		
		term.print("Sorting frames...");
		frames = new Frame[nframes];
		for (int f = 0; f < frames.length; f++) {
			frames[f] = new Frame(f);
		}
		d.reset();
		while (d.hasRemaining()) {
			Trajectory traj = d.nextTrajectory();
			for (int i = 0; i < traj.size(); i++) {
				double t=traj.t(i), x=traj.x(i), y=traj.y(i);
				frames[(int)t].acc(traj.id, x, y);
			}
		}
		term.println("done.");
		
		plot(0);
		Utilities.frame(plot.getComponent(), plot.getTitle());
	}
	
	public FrameAnimator animator() {
		return new FrameAnimator(this);
	}
	
	public Frame getMobileParticles(int t, int tstar, double rstar) {
		int t2 = t;
		int t1 = t2 - tstar;
		if (t1 <= 0) {
			term.println("Too early to print mobile particles");
			return null;
		}
		else {
			Frame f1 = frames[t1];
			Frame f2 = frames[t2];
			Hashtable<Integer,Integer> hash = f1.getIndexHash();
			Frame mobile = new Frame(t2);
			
			for (int i2 = 0; i2 < f2.size(); i2++) {
				int id = (int)f2.id.get(i2);
				
				if (!hash.containsKey(id))
					continue;
				
				int i1 = hash.get(id);
				double x1 = f1.x.get(i1);
				double y1 = f1.y.get(i1);				
				double x2 = f2.x.get(i2);
				double y2 = f2.y.get(i2);				
				
				double d2 = sqr(x2-x1) + sqr(y2-y1);
				if (d2 > sqr(rstar))
					mobile.acc(id, x2, y2);
			}
			return mobile;
		}
	}
	
	boolean bondExists(double x1, double y1, double x2, double y2, double rstar) {
		return (sqr(x2-x1) + sqr(y2-y1) < sqr(rstar));
	}
	
	public List<Frame> getStrings(int t, int tstar, double rstar) {
		ArrayList<Frame> ret = new ArrayList<Frame>();
		
		Frame mobile = getMobileParticles(t, tstar, rstar);
		Hashtable<Integer,Integer> h_pre = frames[t-tstar].getIndexHash();
		
		Percolator<Integer> perc = new Percolator<Integer>();
		for (int i = 0; i < mobile.size(); i++) {
			int id1 = (int)mobile.id.get(i);
			double x1_cur = mobile.x.get(i);
			double y1_cur = mobile.y.get(i);
			double x1_pre = frames[t-tstar].x.get(h_pre.get(id1));
			double y1_pre = frames[t-tstar].y.get(h_pre.get(id1));
			
			perc.add(id1);
			for (int j = 0; j < i; j++) {
				int id2 = (int)mobile.id.get(j);
				double x2_cur = mobile.x.get(j);
				double y2_cur = mobile.y.get(j);
				double x2_pre = frames[t-tstar].x.get(h_pre.get(id2));
				double y2_pre = frames[t-tstar].y.get(h_pre.get(id2));
				
				if (bondExists(x1_cur, y1_cur, x2_pre, y2_pre, rstar) ||
					bondExists(x1_pre, y1_pre, x2_cur, y2_cur, rstar))
					perc.bond(id1, id2);
			}
		}
		
		Hashtable<Integer,Integer> h_cur = mobile.getIndexHash();
		for (ArrayList<Integer> group : perc.getGroups()) {
			Frame f = new Frame(t);
			for (int id : group) {
				f.acc(id, mobile.x.get(h_cur.get(id)), mobile.y.get(h_cur.get(id)));
			}
			ret.add(f);
		}
		return ret;
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
				int particles = frames[t].size();
				double area = bds.getWidth()*bds.getHeight();
				String s = "Frame="+t+" #part="+particles+" len*="+sqrt(area/particles);
				g.drawString(s, 4, 4);
				
				g.setProjection(g.viewBounds());
			}
			public Bounds getBounds() {
				return bds;
			}
		};
	}
	

	public Scene2D plot(int t) {
		if (t >= frames.length)
			term.println("Time " + t + " exceeds maximum frame " + (frames.length-1));
		plot.clearDrawables();
		plot.addDrawable(frames[t].getDrawable(new Color(0f, 0f, 1f, 0.5f)));
		plot.addDrawable(getBackgroundDrawable(t));
		return plot;
	}
	
	public void plotMobile(int t, int tstar, double rstar) {
		if (t >= frames.length)
			term.println("Time " + t + " exceeds maximum frame " + (frames.length-1));
		Frame mobile = getMobileParticles(t, tstar, rstar);
		if (mobile != null) {
			plot.clearDrawables();
			plot.addDrawable(frames[t].getDrawable(new Color(0f, 0f, 1f, 0.3f)));
			plot.addDrawable(mobile.getDrawable(new Color(0f, 0f, 1f, 1f)));
			plot.addDrawable(getBackgroundDrawable(t));
		}
	}
	
	public void plotStrings(int t, List<Frame> strings) {
		if (t >= frames.length)
			term.println("Time " + t + " exceeds maximum frame " + (frames.length-1));

		plot.clearDrawables();
		plot.addDrawable(frames[t].getDrawable(new Color(0f, 0f, 1f, 0.3f)));
		
		Color colors[] = {Color.BLUE, Color.CYAN, Color.GRAY, Color.GREEN,
						  Color.MAGENTA, Color.ORANGE, Color.PINK, Color.RED};
		for (int i = 0; i < strings.size(); i++) {
			Frame f = strings.get(i);
			plot.addDrawable(f.getDrawable(colors[i%colors.length]));
		}
		plot.addDrawable(getBackgroundDrawable(t));
	}
	
	private void calculateSceneBounds(ShakerData d) {
		double tmax = 0;
		bds = new Bounds();
		
		d.reset();
		while (d.hasRemaining()) {
			Trajectory traj = d.nextTrajectory();
			
			for (int i = 0; i < traj.size(); i++) {
				double t=traj.t(i), x=traj.x(i), y=traj.y(i);
				tmax = max(tmax, t);
				bds.xmin = min(bds.xmin, x);
				bds.ymin = min(bds.ymin, y);
				bds.xmax = max(bds.xmax, x);
				bds.ymax = max(bds.ymax, y);
			}
		}
		nframes = 1+(int)tmax;
	}
}

class Frame {
	int t;
	public DynamicArray id = new DynamicArray();
	public DynamicArray x = new DynamicArray();
	public DynamicArray y = new DynamicArray();
	
	public Frame(int t) {
		this.t = t;
	}
	
	public void acc(double id, double x, double y) {
		this.id.append(id);
		this.x.append(x);
		this.y.append(y);
	}
	
	public int size() {
		return id.size();
	}
	
	public Drawable<Gfx2D> getDrawable(final Color color) {
		return new Drawable<Gfx2D>() {
			public void draw(Gfx2D g) {
				g.setColor(color);
				for (int i = 0; i < id.size(); i++) {
					g.fillCircle(x.get(i), y.get(i), 4.0);
				}
			}
			public Bounds getBounds() {
				return new Bounds();
			}
		};
	}
	
	// Maps the particle's global id to the the particle's index in this frame
	public Hashtable<Integer,Integer> getIndexHash() {
		Hashtable<Integer,Integer> ret = new Hashtable<Integer,Integer>();
		for (int i = 0; i < size(); i++) {
			ret.put((int)id.get(i), i);
		}
		return ret;
	}
}


class Trajectory {
	double radius, id;
	public DynamicArray t = new DynamicArray();
	public DynamicArray x = new DynamicArray();
	public DynamicArray y = new DynamicArray();
	
	public void clear() {
		t.clear();
		x.clear();
		y.clear();
	}
	
	public void acc(double x, double y, double t) {
		this.x.append(x);
		this.y.append(y);
		this.t.append(t);
	}
	
	public double dist2(int i1, int i2) {
		double dx = x(i1) - x(i2);
		double dy = y(i1) - y(i2);
		return dx*dx+dy*dy;
	}
	
	public int size() {return t.size();}
	public double x(int i) {return x.get(i); }
	public double y(int i) {return y.get(i); }
	public double t(int i) {return t.get(i); }
}

class ShakerData {
	private int maxCnt = Integer.MAX_VALUE;
	private int cnt;
	private FloatBuffer fb;
	Terminal term;
	
	public ShakerData(Terminal term, String fname) {
		try {
			@SuppressWarnings("resource")
            FileChannel channel = new FileInputStream(fname).getChannel();
			MappedByteBuffer bb = channel.map(FileChannel.MapMode.READ_ONLY, 0, channel.size());
			bb.order(ByteOrder.LITTLE_ENDIAN);
			
			IntBuffer ib = bb.asIntBuffer();
		    ib.get(); // unknown meaning
		    int dim = ib.get(); // spatial dimensions
		    int n = ib.get(); // columns of data
		    int m = ib.get(); // rows of data
		    int prec = ib.get(); // precision (4 for float, 8 for double)
		    int size = ib.get(); // m*n
		    for (int i = 0; i < 6; i++)
		    	ib.get(); // these should be zero
		    assert(dim == 2);
		    assert(prec == 4);
		    assert(size == m*n);
		    
		    bb.position(4*ib.position());
		    
		    cnt = 0;
		    fb = bb.asFloatBuffer();
		    
		} catch (Exception e) {
			term.println(e);
		}
	}
	
	public void reset() {
		cnt = 0;
		fb.rewind();
	}
	
	public boolean hasRemaining() {
		return cnt < maxCnt && fb.hasRemaining();
	}
	
	public Trajectory nextTrajectory() {
		if (!hasRemaining())
			return null;
		
		Trajectory ret = new Trajectory();
		double x, y, radius, time, id;
    	x = fb.get();
    	y = fb.get();
    	fb.get(); // brightness
    	radius = fb.get();
    	time = fb.get();
    	id = fb.get();
    	ret.id = id;
    	ret.radius = radius;
   		ret.acc(x, y, time);
   		
    	while (fb.hasRemaining()) {
    		fb.mark();
        	x = fb.get();
        	y = fb.get();
        	fb.get(); // brightness
        	radius = fb.get();
        	time = fb.get();
        	id = fb.get();
        	if (id == ret.id)
        		ret.acc(x, y, time);
        	else {
        		fb.reset();
        		break;
        	}
    	}
    	cnt++;
		return ret;
	}
}

class Alpha {
	public Accumulator x2;
	public Accumulator x4;
	public Accumulator alpha;
}

class Commands {
	Terminal term;
	int minFrames = 500;
	int frameJump = 100;
	
	public Commands(Terminal term) {
		this.term = term;
	}
	
	public ShakerData loadData(String fname) {
		term.println("Opening " + fname);
		return new ShakerData(term, fname);
	}
	
	public ShakerData loadData() {
		try {
			String fname = FileUtil.loadDialog(term.getConsole(), "Open Shaker Data");
			term.println("Opening " + fname);
			return new ShakerData(term, fname);
		} catch(IOException e) {
			term.println(e);
			return null;
		}
	}
	
	public FrameSequence sequence(ShakerData data) {
		return new FrameSequence(term, data);
	}
	
	@SuppressWarnings("unused")
	public Histogram particlesPerFrame(ShakerData data) {
		int total = 0;
		Histogram len = new Histogram(0.1);
		data.reset();
		while (data.hasRemaining()) {
			Trajectory traj = data.nextTrajectory();
			for (int i = 0; i < traj.size(); i++) {
				len.accum(traj.t(i));
				total += 1;
			}
		}
		return len;
	}
	
	public Histogram trajectoryDistribution(ShakerData data) {
		data.reset();
		Histogram len = new Histogram(100);
		while (data.hasRemaining()) {
			Trajectory traj = data.nextTrajectory();
			len.accum(traj.t.size());
		}
		return len;
	}
	
	public Alpha analyze(ShakerData data) {
		int nsteps = 500;
		int[] steps = new int[nsteps];
		for (int i = 0; i < nsteps; i++) {
			steps[i] = (int)exp(i*0.1);
		}
		
		data.reset();
		double bw = 1;
		Accumulator x2 = new Accumulator(bw);
		Accumulator x4 = new Accumulator(bw);
		while (data.hasRemaining()) {
			Trajectory traj = data.nextTrajectory();
			if (traj.t.size() < minFrames) continue;
			
			for (int i = 0; i < traj.size(); i += frameJump) {
				for (int s = 0; s < nsteps; s++) {
					int j = i+steps[s];
					if (j >= traj.size()) break;
					double d2 = traj.dist2(i, j);
					x2.accum(j-i, d2);
					x4.accum(j-i, d2*d2);
				}
			}
		}
		
		Accumulator alpha = new Accumulator(bw);		
		for (double t : x2.keys()) {
			double d = 2;
			double d2 = x2.eval(t);
			double d4 = x4.eval(t);
			alpha.accum(t, (1/(1.+2./d)) * (d4/(d2*d2)) - 1.0);
		}

		Alpha ret = new Alpha();
		ret.x2 = x2;
		ret.x4 = x4;
		ret.alpha = alpha;
		return ret;
	}
	
	public Histogram vanHove(ShakerData data, int tstar) {
		Histogram ret = new Histogram(0.05);
		ret.setNormalizing(true);
		data.reset();
		while (data.hasRemaining()) {
			Trajectory traj = data.nextTrajectory();
			if (traj.t.size() < minFrames) continue;
			
			for (int i = 0; i < traj.size(); i += frameJump) {
				int j = i + tstar;
				if (j >= traj.size()) break;
				double r = sqrt(traj.dist2(i, j));
				ret.accum(r, 2*PI*r);
			}
		}
		return ret;
	}
	
	public Function vanHoveTheory(int tstar, double diffusion) {
		final double sig2 = diffusion*tstar;
		return new Function(0, 5*sqrt(sig2)) {
			public double eval(double r) {
				return (2*PI*r)*(1/(sig2*2*PI))*exp(-(r*r)/(2*sig2));
			}
		};
	}
}

public class ShakerApp extends Terminal {
	public static void main(String[] args) {
		ShakerApp term = new ShakerApp();
		term.help = "Suggested commands:\n"+
			"(Right click in any plot to change view options or to save data)\n"+
			"\tdata = loadData(); // select tracking data (w/o velocity)\n"+
			"\n"+
			"\t// plot a histogram showing the number of frames over which\n" +
			"\t// individual particles have been tracked\n"+
			"\tplot(trajectoryDistribution(data));\n"+
			"\t// in the following analysis commands, neglect particles tracked\n" +
			"\t// for less than, e.g., 500 frames\n"+
			"\tdata.minFrames = 500; // default 500\n"+
			"\t// smaller frameJump means more averages and slower to process\n"+
			"\tdata.frameJump = 100; // default 100\n"+
			"\ta = analyze(data);\n"+
			"\tplot(a.x2); // coefficient of x2 gives diffusion coef.\n"+
			"\tplot(a.alpha); // maximum of alpha gives t*\n"+
			"\t// do van-hove calculation with t* = 1000\n"+
			"\tplot(vanHove(data, 1000), \"Experiment\");\n"+
			"\t// compare to gaussian with t* = 1000, diffusion coef = 0.003\n"+
			"\t// (the intersection of the two plots gives r*)\n"+
			"\treplot(vanHoveTheory(1000, 0.003), \"Theory\");\n"+
			"\n"+
			"\t// Display the data\n"+
			"\ts = sequence(data); // sequence particles by time\n"+
			"\ts.plot(1); // plot the first frame of the data\n"+
			"\t// plot \"mobile\" particles using t=1000, t*=100, r*=10\n"+
			"\ts.plotMobile(1000, 100, 10);\n"+
			"\n"+
			"\t// String analysis\n"+
			"\t// Get a list of strings using t = 1000, t*=100, r*=10\n"+
			"\tstrs = s.getStrings(1000, 100, 10);\n"+
			"\ts.plotStrings(1000, strs); // plot strings at t=1000\n"+
			"\n"+
			"\t// Animate the data\n"+
			"\ta = s.animator();\n" +
			"\ta.start(); // animate the tracked particles\n"+
			"\ta.stop();\n" +
			"\ta.setMobileParams(100, 10); // set (t*, r*) for mobile particles\n"+
			"\ta.selectImageDir(); // (optional) write animation frames to images\n"+
			"\ta.start(); // animate mobile particles\n";
		term.importObject(new Commands(term));
		term.runApplication();
	}
}
